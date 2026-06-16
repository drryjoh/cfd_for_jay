#include "codeJeNN_diffusion.H"

// ============================================================
// nn_diffusion_interface.H
//
// Interface between OpenFOAM and the CodeJeNN-generated
// diffusion coefficient NN (codeJeNN_diffusion.H).
//
// Mirrors the NNInterface pattern from the chemistry NN
// (old_codejenn/nn_interface.H) but targets mixture-averaged
// diffusion coefficients instead of reaction source terms.
// ============================================================

struct NNDiffInterface
{
    const fvMesh& mesh_;

    // Constructor: stores the mesh reference.
    // Call check_species_list(Y) separately during initialisation
    // when useDiffNN is true (mirrors the NNInterface pattern).
    NNDiffInterface(const fvMesh& mesh)
    :
        mesh_(mesh)
    {}

    // Return a list of the species names the NN expects,
    // in the order the NN was trained on.
    static Foam::stringList species_names()
    {
        Foam::stringList names(n_species_diff);
        forAll(names, i)
            names[i] = species_name_diff_gen(i);
        return names;
    }

    static Foam::string species_name(const Foam::label specieI)
    {
        return species_name_diff_gen(specieI);
    }

    // Verify that Y[i].name() matches the NN's expected species
    // order.  Called once during solver initialisation so any
    // mismatch is caught before the time loop begins.
    void check_species_list(const PtrList<volScalarField>& Y) const
    {
        if (Y.size() != n_species_diff)
        {
            FatalErrorInFunction
                << "Species count mismatch: solver has " << Y.size()
                << " species, NN expects " << n_species_diff << "."
                << exit(FatalError);
        }

        Info << "NNDiffInterface: checking species order" << endl;
        Foam::stringList names = species_names();
        forAll(names, i)
        {
            if (Y[i].name() != names[i])
            {
                FatalErrorInFunction
                    << "Species mismatch at index " << i
                    << ": solver has '" << Y[i].name()
                    << "', NN expects '" << names[i] << "'."
                    << exit(FatalError);
            }
        }
        Info << "NNDiffInterface: species check passed" << endl;
    }

    // Populate DiffNN[i] (m^2/s) for every cell and species
    // using the NN given the current T, p, and mass fractions Y.
    void compute(
        PtrList<volScalarField>& DiffNN,
        const PtrList<volScalarField>& Y,
        const volScalarField& p,
        const volScalarField& T) const
    {
        static const scalar inv_patm = scalar(1) / scalar(101325);

        std::array<scalar, n_species_diff + 2> nn_input;

        forAll(T, celli)
        {
            nn_input[0] = T[celli];
            nn_input[1] = p[celli] * inv_patm;
            for (label i = 0; i < n_species_diff; ++i)
                nn_input[2 + i] = Y[i][celli];

            const auto output = predict_diffusion(nn_input);

            for (label i = 0; i < n_species_diff; ++i)
                DiffNN[i][celli] = Foam::max(output[i], scalar(SMALL));
        }
    }
};
