/**
 * Construction d'une replique (reseau + maille) a partir d'un fichier {@link Input}.
 *
 * <p>La recette — ordre des sites de base, produit {@code transpose(K) . g . K} qui donne les
 * g-tenseurs, {@code setGtensor(g0)}, le sens des interactions {@code (b1, b2)} pour le
 * sous-reseau 0 et {@code (b2, b1)} pour le sous-reseau 1, puis {@code createTemplates} avant
 * {@code setInteractionSites} — etait historiquement recopiee a l'identique dans
 * {@code AppParallel.replica} et dans {@code MonteCarloParallelTest.buildFromInput}. Trois copies,
 * dont une dans un test dont le role est precisement de verifier que le chemin parallele construit
 * le <i>meme</i> hamiltonien que la production : si la copie du test derive, le test continue de
 * passer et ne prouve plus rien. La recette vit donc ici, et le test la traverse comme la
 * production.</p>
 *
 * <h2>Champ magnetique</h2>
 * <p>{@link #build(Input, int, int)} installe un champ <b>nul</b>. Le champ depend de l'amplitude
 * h de la boucle et de la direction choisie ; il se calcule a partir du g-tenseur du reseau, qui
 * n'existe qu'une fois la construction faite. L'appelant fait donc :</p>
 * <pre>
 *   ReplicaFactory.Built b = ReplicaFactory.build(in, La, Lb);
 *   double[] H = MathOps.scale(h, MathOps.matrixXvector(b.lattice().getGtensor(),
 *                                                     MathOps.normalize(dir)));
 *   b.lattice().setInteractionField(H);
 * </pre>
 * <p>{@code setInteractionField} ne fait que stocker la reference : elle peut etre appelee a
 * tout moment, avant comme apres {@code createTemplates} (les gabarits ne contiennent que
 * l'echange).</p>
 */
public final class ReplicaFactory {

    private ReplicaFactory() {}

    /**
     * Reseau et maille d'une replique. Les deux sont necessaires : {@code uc} est reclame par
     * {@code DynamicStructureFactor.processSQW}, {@code lattice} par tout le reste.
     */
    public record Built(Lattice lattice, UnitCell uc) {}

    /**
     * Construit la replique.
     *
     * @param in fichier d'entree deja lu.
     * @param La nombre de mailles suivant a1.
     * @param Lb nombre de mailles suivant a2.
     * @return le reseau {@code La x Lb x 2} et sa maille, champ magnetique nul.
     * @throws IllegalArgumentException si {@code in} est {@code null} ou si {@code La} ou
     *         {@code Lb} est &lt; 1.
     */
    public static Built build(Input in, int La, int Lb) {
        if (in == null) throw new IllegalArgumentException("in == null");
        if (La < 1 || Lb < 1) {
            throw new IllegalArgumentException("La et Lb doivent etre >= 1 : " + La + ", " + Lb);
        }

        // Vecteurs primitifs et maille (memes conventions que App et HoneycombGeometry).
        final UnitCell uc = new UnitCell(HoneycombGeometry.a1(), HoneycombGeometry.a2());

        // Positions des deux sites de base ; addBasisSite renvoie leur index.
        final int b1 = uc.addBasisSite(in.getAtomicPositions()[0]);
        final int b2 = uc.addBasisSite(in.getAtomicPositions()[1]);

        final Lattice lattice = new Lattice(new int[] { La, Lb, HoneycombGeometry.SUBLATTICES });

        // Bases de Kitaev a*bc et abc.
        lattice.setKitaevBasis(in.getKitaevBasis()[0]);
        lattice.setCOB(in.getChangeOfBasis());

        // g-tenseurs exprimes dans la base de Kitaev : transpose(K) . g . K.
        final double[][] g0 = MathOps.mult2D(MathOps.transpose(in.getKitaevBasis()[0]),
                MathOps.mult2D(in.getGTensors()[0], in.getKitaevBasis()[0]));
        final double[][] g1 = MathOps.mult2D(MathOps.transpose(in.getKitaevBasis()[1]),
                MathOps.mult2D(in.getGTensors()[1], in.getKitaevBasis()[1]));
        uc.setInteractionGtensor(b1, g0);
        uc.setInteractionGtensor(b2, g1);
        lattice.setGtensor(g0);

        // Champ nul par defaut ; l'appelant l'installe apres coup (cf. javadoc de classe).
        lattice.setInteractionField(new double[3]);

        // Couplages spin-spin : sous-reseau 0 vers 1, puis 1 vers 0.
        for (int sr = 0; sr < lattice.lengthSr; sr++) {
            for (int i = 0; i < in.getInteractionSites()[sr].length; i++) {
                if (sr == 0) {
                    uc.addInteraction(b1, b2, in.getInteractionMatrices()[0][i],
                            in.getInteractionSites()[0][i]);
                }
                if (sr == 1) {
                    uc.addInteraction(b2, b1, in.getInteractionMatrices()[1][i],
                            in.getInteractionSites()[1][i]);
                }
            }
        }

        lattice.createTemplates(uc);
        lattice.setInteractionSites(uc);

        // Terme quartique optionnel (premiers voisins) lu depuis le fichier d'entree ;
        // absent (b = 0) -> aucune physique ajoutee, reseau identique a l'historique.
        if (in.getQuarticB() != 0.d) {
            lattice.setQuartic(in.getQuarticB(), in.getQuarticGammaSub0(), in.getQuarticGammaSub1());
        }

        return new Built(lattice, uc);
    }
}
