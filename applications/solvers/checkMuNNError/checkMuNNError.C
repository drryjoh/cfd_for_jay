/*---------------------------------------------------------------------------*\
  checkMuNNError.C

  Post-processing utility for detonationFoam NS_mixtureAverage cases.

  Reads the solver-computed mixMu field and species/temperature fields from
  the latest time directory, evaluates NNHI and NNLO, and writes the
  cell-by-cell relative error for each model:

      mu_err_hi   |mixMu - mu_NNHI| / mixMu
      mu_err_lo   |mixMu - mu_NNLO| / mixMu

  Using mixMu directly (written by the solver via AUTO_WRITE) means the
  reference is always the exact Wilke value the solver used, with no risk
  of coefficient or unit mismatch.

  Fields read from the latest time directory:
      mixMu, H2, O2, N2, T

  NN input order: [Y_O2, Y_N2, Y_H2, T]   (matches nn_mu_interface.H)
  Species order:  Y[0]=H2, Y[1]=O2, Y[2]=N2

  Usage (run from the case directory):
      checkMuNNError
\*---------------------------------------------------------------------------*/

#include "model_hi_NEW.hpp"   // NNHI: model_hi_NEW([Y_O2, Y_N2, Y_H2, T])
#include "model_lo_NEW.hpp"   // NNLO: model_lo_NEW([Y_O2, Y_N2, Y_H2, T])

#include "fvCFD.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();

    #include "setRootCaseLists.H"
    #include "createTime.H"

    // Set runTime to the latest available time directory
    instantList times = runTime.times();
    if (times.empty())
        FatalErrorInFunction
            << "No time directories found in " << runTime.path()
            << exit(FatalError);

    runTime.setTime(times.last(), times.size() - 1);
    Info << "Processing time = " << runTime.timeName() << nl << endl;

    #include "createMesh.H"

    // ----------------------------------------------------------------
    // Read fields from the latest time directory
    // ----------------------------------------------------------------
    volScalarField mixMu
    (
        IOobject("mixMu", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );
    volScalarField Y_H2
    (
        IOobject("H2",    runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );
    volScalarField Y_O2
    (
        IOobject("O2",    runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );
    volScalarField Y_N2
    (
        IOobject("N2",    runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );
    volScalarField T
    (
        IOobject("T",     runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );

    // ----------------------------------------------------------------
    // Output error fields
    // ----------------------------------------------------------------
    volScalarField mu_err_hi
    (
        IOobject
        (
            "mu_err_hi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, scalar(0))
    );

    volScalarField mu_err_lo
    (
        IOobject
        (
            "mu_err_lo",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, scalar(0))
    );

    // ----------------------------------------------------------------
    // Cell-by-cell computation
    // ----------------------------------------------------------------
    scalar sumErrHI = 0, sumErrLO = 0;
    scalar maxErrHI = 0, maxErrLO = 0;

    forAll(T, celli)
    {
        // NN input order: [O2, N2, H2, T]
        std::array<scalar, 4> nn_in =
        {
            Y_O2[celli], Y_N2[celli], Y_H2[celli], T[celli]
        };

        const scalar muRef = mixMu[celli];
        const scalar muHI  = std::max(model_hi_NEW(nn_in)[0], scalar(1e-30));
        const scalar muLO  = std::max(model_lo_NEW(nn_in)[0], scalar(1e-30));

        const scalar eHI = std::abs(muRef - muHI) / muRef;
        const scalar eLO = std::abs(muRef - muLO) / muRef;

        mu_err_hi[celli] = eHI;
        mu_err_lo[celli] = eLO;

        sumErrHI += eHI;  maxErrHI = std::max(maxErrHI, eHI);
        sumErrLO += eLO;  maxErrLO = std::max(maxErrLO, eLO);
    }

    mu_err_hi.write();
    mu_err_lo.write();

    const label nCells = mesh.nCells();
    Info << "Written to " << runTime.timeName() << nl
         << "  NNHI | mean: " << sumErrHI/nCells*100 << " %"
                   << "  max: " << maxErrHI*100 << " %" << nl
         << "  NNLO | mean: " << sumErrLO/nCells*100 << " %"
                   << "  max: " << maxErrLO*100 << " %" << endl;

    return 0;
}
