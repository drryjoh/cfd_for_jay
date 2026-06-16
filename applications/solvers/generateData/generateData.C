/*---------------------------------------------------------------------------*\
  generateData.C

  Standalone data-generation utility for training viscosity NNs.

  For each sample, one of two sampling modes is used:
    a) Uniform Dirichlet  (fraction: 1 - fracPhysical)
         Y_i ~ Uniform(0,1), normalised to sum = 1.
         T   ~ Uniform(200, 3000) K.
    b) Physical mixing-line  (fraction: fracPhysical)
         Y = alpha * Y_fuel + (1-alpha) * Y_oxidizer,  alpha ~ Uniform(0,1).
         T = alpha * T_fuel + (1-alpha) * T_oxidizer  +  scatter.
         This directly covers the inlet-to-inlet composition trajectory that
         appears in splitter-plate / mixing-layer simulations and is
         underrepresented by the Dirichlet distribution (which concentrates
         near equal-composition mixtures).

  For each state:
    1. Compute per-species viscosity from log-polynomial fits.
    2. Apply the Wilke mixture rule (identical to updateTransProperties.H).
    3. Write one CSV line:  j, Y_0, ..., Y_{n-1}, T, mu_mix

  Run from a case directory that contains:
      constant/generateDataProperties

  Usage:
      generateData
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"

#include <fstream>
#include <random>
#include <cmath>

using namespace Foam;

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
        props.lookupOrDefault<word>("outputFile", "mu_training_data.csv");

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
        W[i]    = readScalar(sd.lookup("W"));
        mu1[i]  = readScalar(sd.lookup("Mu1"));
        mu2[i]  = readScalar(sd.lookup("Mu2"));
        mu3[i]  = readScalar(sd.lookup("Mu3"));
        mu4[i]  = readScalar(sd.lookup("Mu4"));
    }

    // ----------------------------------------------------------------
    // Physical mixing-line importance sampling
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
            Y_fuel[0] = 1.0;   // default: pure first species
        }

        if (props.found("Y_oxidizer"))
        {
            List<scalar> tmp(props.lookup("Y_oxidizer"));
            Y_ox = tmp;
        }
        else if (nSp >= 3)
        {
            Y_ox[1] = 0.233;   // default: air-like O2
            Y_ox[2] = 0.767;   // default: air-like N2
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

    Info << "Generating " << nSamples << " samples"
         << " | species: " << speciesNames << nl << endl;

    // ----------------------------------------------------------------
    // Random sampling
    // ----------------------------------------------------------------
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<scalar> uni(0.0, 1.0);

    std::ofstream out(outputFile);

    // Header
    out << "j";
    forAll(speciesNames, i)
        out << "," << speciesNames[i];
    out << ",T,mu\n";

    for (label j = 0; j < nSamples; ++j)
    {
        List<scalar> Y(nSp);
        scalar T;

        if (fracPhysical > 0 && uni(rng) < fracPhysical)
        {
            // -- Physical mixing-line sample ---------------------------
            // alpha=1 → pure fuel;  alpha=0 → pure oxidizer
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

        // -- Per-species viscosity [Pa·s] --------------------------
        List<scalar> Mu(nSp);
        for (int i = 0; i < nSp; ++i)
            Mu[i] = 0.1 * std::exp(
                mu1[i] + logT*(mu2[i] + logT*(mu3[i] + logT*mu4[i])));

        // -- Mole fractions ----------------------------------------
        scalar invWmix = 0;
        for (int i = 0; i < nSp; ++i) invWmix += Y[i] / W[i];
        const scalar Wmix = 1.0 / invWmix;

        List<scalar> X(nSp);
        for (int i = 0; i < nSp; ++i) X[i] = Y[i] * Wmix / W[i];

        // -- Wilke mixture rule ------------------------------------
        static const scalar inv_sqrt8 = 1.0 / std::sqrt(8.0);
        scalar mixMu = 0;
        for (int i = 0; i < nSp; ++i)
        {
            scalar denom = 0;
            for (int z = 0; z < nSp; ++z)
            {
                denom += X[z] * inv_sqrt8
                       * std::pow(1.0 + std::sqrt(Mu[i]/Mu[z])
                                      * std::pow(W[z]/W[i], 0.25), 2.0)
                       / std::sqrt(1.0 + W[i]/W[z]);
            }
            mixMu += X[i] * Mu[i] / denom;
        }

        // -- Write one line ----------------------------------------
        out << j;
        for (int i = 0; i < nSp; ++i)
            out << "," << Y[i];
        out << "," << T << "," << mixMu << "\n";
    }

    out.close();

    Info << "Done.  " << nSamples << " samples written to '"
         << outputFile << "'" << endl;

    return 0;
}
