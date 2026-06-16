/*---------------------------------------------------------------------------*\
  checkNNPrediction.C

  Validates two CodeJeNN viscosity NNs (NNHI, NNLO) against the
  polynomial Wilke model.

  Sampling matches generateData.C exactly:
    a) Uniform Dirichlet  (fraction: 1 - fracPhysical)
    b) Physical mixing-line  (fraction: fracPhysical)
         Y = alpha * Y_fuel + (1-alpha) * Y_oxidizer,  alpha ~ Uniform(0,1)
         T blended + scatter

  For each state:
    1. Compute polynomial per-species Mu_i, apply Wilke rule → mu_wilke.
    2. Call each NN model → mu_NNHI, mu_NNLO.
    3. Relative error = |mu_wilke - mu_NNk| / mu_wilke.
    4. Write CSV + print summary statistics.

  Reads from:  constant/generateDataProperties
  Writes to:   outputFile (default: nn_check_results.csv)

  Usage:
      checkNNPrediction
\*---------------------------------------------------------------------------*/

#include "model_hi_NEW.hpp"   // NNHI: model_hi_NEW([Y_O2, Y_N2, Y_H2, T])
#include "model_lo_NEW.hpp"   // NNLO: model_lo_NEW([Y_O2, Y_N2, Y_H2, T])

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"

#include <fstream>
#include <random>
#include <cmath>
#include <limits>

using namespace Foam;

// Wilke mixture viscosity rule.
static scalar wilkeMix
(
    const List<scalar>& Mu,
    const List<scalar>& X,
    const List<scalar>& W
)
{
    static const scalar inv_sqrt8 = 1.0 / std::sqrt(8.0);
    const int n = Mu.size();
    scalar mix = 0;
    for (int i = 0; i < n; ++i)
    {
        scalar denom = 0;
        for (int z = 0; z < n; ++z)
        {
            denom += X[z] * inv_sqrt8
                   * std::pow(1.0 + std::sqrt(Mu[i]/Mu[z])
                                  * std::pow(W[z]/W[i], 0.25), 2.0)
                   / std::sqrt(1.0 + W[i]/W[z]);
        }
        mix += X[i] * Mu[i] / denom;
    }
    return mix;
}


int main(int argc, char *argv[])
{
    argList::noParallel();

    #include "setRootCaseLists.H"
    #include "createTime.H"

    // ----------------------------------------------------------------
    // Read configuration
    // ----------------------------------------------------------------
    IOdictionary props
    (
        IOobject
        (
            "generateDataProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const label nSamples =
        readLabel(props.lookup("nSamples"));

    const word outputFile =
        props.lookupOrDefault<word>("outputFile", "nn_check_results.csv");

    // ----------------------------------------------------------------
    // Species data
    // ----------------------------------------------------------------
    const wordList speciesNames(props.lookup("species"));
    const int nSp = speciesNames.size();

    List<scalar> W(nSp);
    List<scalar> mu1(nSp), mu2(nSp), mu3(nSp), mu4(nSp);

    forAll(speciesNames, i)
    {
        const dictionary& sd = props.subDict(speciesNames[i]);
        W[i]   = readScalar(sd.lookup("W"));
        mu1[i] = readScalar(sd.lookup("Mu1"));
        mu2[i] = readScalar(sd.lookup("Mu2"));
        mu3[i] = readScalar(sd.lookup("Mu3"));
        mu4[i] = readScalar(sd.lookup("Mu4"));
    }

    // ----------------------------------------------------------------
    // Physical mixing-line importance sampling (matches generateData.C)
    // ----------------------------------------------------------------
    const scalar fracPhysical =
        props.lookupOrDefault<scalar>("fracPhysical", 0.0);

    List<scalar> Y_fuel(nSp, scalar(0));
    List<scalar> Y_ox(nSp, scalar(0));
    scalar T_fuel    = 350.0;
    scalar T_ox      = 1500.0;
    scalar T_scatter = 200.0;

    if (fracPhysical > 0)
    {
        if (props.found("Y_fuel"))
        {
            List<scalar> tmp(props.lookup("Y_fuel"));
            Y_fuel = tmp;
        }
        else
        {
            Y_fuel[0] = 1.0;
        }

        if (props.found("Y_oxidizer"))
        {
            List<scalar> tmp(props.lookup("Y_oxidizer"));
            Y_ox = tmp;
        }
        else if (nSp >= 3)
        {
            Y_ox[1] = 0.233;
            Y_ox[2] = 0.767;
        }

        T_fuel    = props.lookupOrDefault<scalar>("T_fuel",      350.0);
        T_ox      = props.lookupOrDefault<scalar>("T_oxidizer", 1500.0);
        T_scatter = props.lookupOrDefault<scalar>("T_scatter",   200.0);

        Info << "Physical mixing-line sampling enabled: fracPhysical = "
             << fracPhysical << nl
             << "  Y_fuel    = " << Y_fuel    << nl
             << "  Y_oxidizer= " << Y_ox      << nl
             << "  T_fuel    = " << T_fuel    << " K" << nl
             << "  T_oxidizer= " << T_ox      << " K" << nl
             << "  T_scatter = " << T_scatter << " K" << endl;
    }

    Info << "Checking NN predictions for " << nSamples << " samples"
         << " | species: " << speciesNames << nl << endl;

    // ----------------------------------------------------------------
    // Random sampling
    // ----------------------------------------------------------------
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<scalar> uni(0.0, 1.0);

    std::ofstream out(outputFile);

    out << "j,Y_H2,Y_O2,Y_N2,T,mu_wilke,mu_NNHI,mu_NNLO,eHI,eLO\n";

    scalar sumErrHI = 0, sumErrLO = 0;
    scalar maxErrHI = 0, maxErrLO = 0;

    for (label j = 0; j < nSamples; ++j)
    {
        List<scalar> Y(nSp);
        scalar T;

        if (fracPhysical > 0 && uni(rng) < fracPhysical)
        {
            // -- Physical mixing-line sample ---------------------------
            const scalar alpha = uni(rng);
            for (int i = 0; i < nSp; ++i)
                Y[i] = alpha * Y_fuel[i] + (1.0 - alpha) * Y_ox[i];

            const scalar T_blend = alpha * T_fuel + (1.0 - alpha) * T_ox;
            T = T_blend + (2.0 * uni(rng) - 1.0) * T_scatter;
            T = std::max(T, scalar(200.0));
            T = std::min(T, scalar(3000.0));
        }
        else
        {
            // -- Uniform Dirichlet + full T range ----------------------
            scalar Ysum = 0;
            for (int i = 0; i < nSp; ++i) { Y[i] = uni(rng); Ysum += Y[i]; }
            for (int i = 0; i < nSp; ++i) Y[i] /= Ysum;
            T = 200.0 + uni(rng) * 2800.0;
        }

        const scalar logT = std::log(T);

        // -- Mole fractions ----------------------------------------
        scalar invWmix = 0;
        for (int i = 0; i < nSp; ++i) invWmix += Y[i] / W[i];
        const scalar Wmix = 1.0 / invWmix;

        List<scalar> X(nSp);
        for (int i = 0; i < nSp; ++i) X[i] = Y[i] * Wmix / W[i];

        // -- Polynomial per-species Mu [Pa·s] ----------------------
        List<scalar> MuPoly(nSp);
        for (int i = 0; i < nSp; ++i)
            MuPoly[i] = 0.1 * std::exp(
                mu1[i] + logT*(mu2[i] + logT*(mu3[i] + logT*mu4[i])));

        const scalar muWilke = wilkeMix(MuPoly, X, W);

        // -- NN predictions [Pa·s] ---------------------------------
        // Species order in generateDataProperties: H2=Y[0], O2=Y[1], N2=Y[2]
        // NN input order: [O2, N2, H2, T]
        std::array<scalar, 4> nn_input = { Y[1], Y[2], Y[0], T };

        const scalar muHI = std::max(model_hi_NEW(nn_input)[0], scalar(1e-30));
        const scalar muLO = std::max(model_lo_NEW(nn_input)[0], scalar(1e-30));

        // -- Relative errors ---------------------------------------
        const scalar eHI = std::abs(muWilke - muHI) / muWilke;
        const scalar eLO = std::abs(muWilke - muLO) / muWilke;

        sumErrHI += eHI;  maxErrHI = std::max(maxErrHI, eHI);
        sumErrLO += eLO;  maxErrLO = std::max(maxErrLO, eLO);

        // -- Write -------------------------------------------------
        out << j
            << "," << Y[0] << "," << Y[1] << "," << Y[2]
            << "," << T
            << "," << muWilke
            << "," << muHI << "," << muLO
            << "," << eHI  << "," << eLO
            << "\n";
    }

    out.close();

    Info << "Results written to '" << outputFile << "'" << nl
         << "  NNHI | mean: " << sumErrHI/nSamples*100 << " %"
                   << "  max: " << maxErrHI*100 << " %" << nl
         << "  NNLO | mean: " << sumErrLO/nSamples*100 << " %"
                   << "  max: " << maxErrLO*100 << " %" << endl;

    return 0;
}
