#include "model_hi_NEW.hpp"
#include "model_lo_NEW.hpp"

// ============================================================
// nn_mu_interface.H
//
// Interface between OpenFOAM and the CodeJeNN viscosity NNs.
// Two models are selectable via muNNmodel in chemistryProperties:
//
//   NNHI — model_hi_NEW.hpp  (high-accuracy, deeper architecture)
//   NNLO — model_lo_NEW.hpp  (low-cost, shallow architecture)
//
// Input order: [Y_O2, Y_N2, Y_H2, T]  (matches training: df[['O2','N2','H2','T']])
// Solver species order: Y[0]=H2, Y[1]=O2, Y[2]=N2.
// ============================================================

// Species count and name table used by check_species_list.
const int n_species_mu = 3;

static std::string species_name_mu_gen(int index)
{
    constexpr std::array<const char*, 3> names = {"H2", "O2", "N2"};
    return names[index];
}

struct NNMuInterface
{
    const fvMesh& mesh_;
    label modelIdx_;

    NNMuInterface(const fvMesh& mesh, const word& modelName)
    :
        mesh_(mesh),
        modelIdx_(nameToIdx(modelName))
    {
        Info << "NNMuInterface: active model = " << modelName << endl;
    }

    static label nameToIdx(const word& name)
    {
        if (name == "NNHI") return 0;
        if (name == "NNLO") return 1;
        FatalErrorInFunction
            << "Unknown muNNmodel '" << name
            << "'. Valid options: NNHI, NNLO."
            << exit(FatalError);
        return 0;
    }

    void check_species_list(const PtrList<volScalarField>& Y) const
    {
        if (Y.size() != n_species_mu)
        {
            FatalErrorInFunction
                << "Species count mismatch: solver has " << Y.size()
                << " species, Mu NN expects " << n_species_mu << "."
                << exit(FatalError);
        }

        Info << "NNMuInterface: checking species order" << endl;
        for (int i = 0; i < n_species_mu; ++i)
        {
            if (Y[i].name() != species_name_mu_gen(i))
            {
                FatalErrorInFunction
                    << "Species mismatch at index " << i
                    << ": solver has '" << Y[i].name()
                    << "', Mu NN expects '" << species_name_mu_gen(i) << "'."
                    << exit(FatalError);
            }
        }
        Info << "NNMuInterface: species check passed" << endl;
    }

    // Raw NN output [Pa·s] for a single cell — no clamping.
    // Returns a negative value when the NN extrapolates out of range.
    inline scalar predictCellRaw
    (
        scalar T_c,
        scalar Y_H2,
        scalar Y_O2,
        scalar Y_N2
    ) const
    {
        // Input order: [O2, N2, H2, T]
        std::array<scalar, 4> in = {Y_O2, Y_N2, Y_H2, T_c};
        switch (modelIdx_)
        {
            case 0: return model_hi_NEW(in)[0];
            case 1: return model_lo_NEW(in)[0];
            default: return scalar(0);
        }
    }

    // Clamped prediction — returns SMALL when raw output is non-positive.
    inline scalar predictCell
    (
        scalar T_c,
        scalar Y_H2,
        scalar Y_O2,
        scalar Y_N2
    ) const
    {
        return Foam::max(predictCellRaw(T_c, Y_H2, Y_O2, Y_N2), scalar(SMALL));
    }

    // Fill mixMu [Pa·s] for every cell using the active NN.
    void compute
    (
        volScalarField& mixMu,
        const PtrList<volScalarField>& Y,
        const volScalarField& T
    ) const
    {
        forAll(T, celli)
        {
            mixMu[celli] = predictCell
            (
                T[celli],
                Y[0][celli],   // H2
                Y[1][celli],   // O2
                Y[2][celli]    // N2
            );
        }
    }

    // Fill mixMu and increment errorCount for cells where the raw NN output
    // was non-positive (i.e., physically invalid — clamped to SMALL).
    void compute
    (
        volScalarField& mixMu,
        const PtrList<volScalarField>& Y,
        const volScalarField& T,
        volScalarField& errorCount
    ) const
    {
        forAll(T, celli)
        {
            const scalar rawMu = predictCellRaw
            (
                T[celli],
                Y[0][celli],   // H2
                Y[1][celli],   // O2
                Y[2][celli]    // N2
            );

            if (rawMu < scalar(SMALL))
            {
                mixMu[celli] = scalar(SMALL);
                errorCount[celli] += scalar(1);
            }
            else
            {
                mixMu[celli] = rawMu;
            }
        }
    }
};
