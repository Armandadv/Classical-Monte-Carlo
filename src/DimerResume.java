import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.avro.file.DataFileReader;
import org.apache.avro.generic.GenericDatumReader;
import org.apache.avro.generic.GenericRecord;

/**
 * Reprise du Monte Carlo a la derniere temperature + recalcul du S_dimer connecte,
 * run par run, SANS aucune autre sortie.
 *
 * <p>Pour CHAQUE dossier de run sous la racine (un dossier = contenant des
 * {@code spins_T*.avro}) l'outil :
 * <ol>
 *   <li>selectionne le groupe d'instantanes a la plus basse temperature (a 1e-3 pres) ;</li>
 *   <li>retrouve le fichier d'input du run dans le dossier des inputs, par nom : le
 *       dossier {@code 13092026KBQ_1 - Copie (278)_h0.000000} correspond au fichier
 *       {@code KBQ_1 - Copie (278)} (date preferee et suffixe {@code _h<valeur>} retires) ;</li>
 *   <li>construit le reseau (RecetteFactory/Input), y charge chaque instantane, relance
 *       le MC a T fixe (thermalisation puis nMesures mesures separees de rate sweeps,
 *       sigma adaptatif comme en production) et accumule S_dimer ;</li>
 *   <li>ecrit UNIQUEMENT {@code structure_dimere_<nomRun>_T<temp>_h0.000000.avro} dans --out. Aucune ecriture dans les dossiers de run : pas
 *       d'equilibration CSV, pas d'observables, pas de SQW — les merges existants
 *       ne voient rien de nouveau.</li>
 * </ol></p>
 *
 * <p>Convention S_dimer : correlateur CONNECTE de {@link DimerStructureFactor} (fond
 * retire par configuration), cf. Eq. (2) de arXiv:2202.11100 (variante fond total).
 * Champ magnetique : seuls les runs {@code _h0.000000} sont traites (le champ n'est
 * pas installe par cet outil — les dossiers h != 0 sont signales et sautes).</p>
 *
 * <h2>Usage</h2>
 * <pre>
 * java -cp "bin:lib/*" DimerResume &lt;racine_runs&gt; &lt;dossier_inputs&gt; &lt;La&gt; &lt;Lb&gt;
 *        [--h 3] [--k 3] [--nmes 200] [--rate 20] [--therm n] [--sig 0.19]
 *        [--overrelax 0] [--threads 1] [--out .] [--seed n] [--temp T]
 * </pre>
 * {@code --temp T} : au lieu du groupe le plus froid, reprend le MC au groupe de
 * temperature le plus proche de T (chasse au pinch point dans la fenetre liquide,
 * au-dessus de la transition).
 */
public class DimerResume {

    private static final double TOL_T = 1e-3;
    private static final Pattern TAG_T = Pattern.compile("^spins_T([-+0-9.eE]+)_");
    private static final Pattern SUFFIXE_H = Pattern.compile("_h([-+0-9.eE]+)$");

    public static void main(String[] args) throws IOException {
        String racine = null, inputsDir = null, outDir = ".";
        int La = 0, Lb = 0, H = 3, K = 3, nMes = 200, rate = 20, therm = -1;
        int overrelax = 0, threads = 1;
        double sig0 = 0.19;
        double tempCible = Double.NaN;
        long seed = System.nanoTime();

        for (int i = 0; i < args.length; i++) {
            final String a = args[i];
            switch (a) {
                case "--h": H = Integer.parseInt(args[++i]); break;
                case "--k": K = Integer.parseInt(args[++i]); break;
                case "--nmes": nMes = Integer.parseInt(args[++i]); break;
                case "--rate": rate = Integer.parseInt(args[++i]); break;
                case "--therm": therm = Integer.parseInt(args[++i]); break;
                case "--sig": sig0 = Double.parseDouble(args[++i]); break;
                case "--overrelax": overrelax = Integer.parseInt(args[++i]); break;
                case "--threads": threads = Integer.parseInt(args[++i]); break;
                case "--out": outDir = args[++i]; break;
                case "--seed": seed = Long.parseLong(args[++i]); break;
                case "--temp": tempCible = Double.parseDouble(args[++i]); break;
                default:
                    if (racine == null) racine = a;
                    else if (inputsDir == null) inputsDir = a;
                    else if (La == 0) La = Integer.parseInt(a);
                    else if (Lb == 0) Lb = Integer.parseInt(a);
                    else throw new IllegalArgumentException("argument en trop : " + a);
            }
        }
        if (racine == null || inputsDir == null || La == 0 || Lb == 0) {
            System.err.println("Usage : java DimerResume <racine_runs> <dossier_inputs> "
                    + "<La> <Lb> [--h 3] [--k 3] [--nmes 200] [--rate 20] [--therm n] "
                    + "[--sig 0.19] [--overrelax 0] [--threads 1] [--out .] [--seed n]");
            System.exit(1);
        }
        if (therm < 0) therm = rate;

        final File root = new File(racine);
        final List<File> dossiers = dossiersDeRun(root);
        System.out.println(dossiers.size() + " dossier(s) de run sous " + racine);

        int ok = 0, sautes = 0;
        for (final File dossier : dossiers) {
            try {
                traiterDossier(dossier, new File(inputsDir), La, Lb, H, K,
                        nMes, rate, therm, overrelax, threads, sig0, seed, outDir,
                        tempCible);
                ok++;
            } catch (IOException | RuntimeException e) {
                // RuntimeException aussi : Input(String) avale les IOException de lecture
                // (demi-objet construit sur un fichier malforme) et ReplicaFactory.build
                // leverait alors un NPE — un dossier pourri ne doit pas tuer le batch.
                sautes++;
                System.out.println("  SKIP " + dossier.getName() + " : " + e.getMessage());
            }
        }
        System.out.println("\ntermine : " + ok + " dossier(s) traites, " + sautes
                + " saute(s). Sorties dans " + new File(outDir).getAbsolutePath());
    }

    /** Erreur non fatale : le dossier est saute avec un message. */
    private static class SauteException extends IOException {
        SauteException(String msg) { super(msg); }
    }

    // ------------------------------------------------------------------

    private static List<File> dossiersDeRun(File root) {
        final List<File> res = new ArrayList<>();
        final File[] subs = root.listFiles(File::isDirectory);
        if (subs != null) {
            for (final File d : subs) {
                final String[] l = d.list((dir, nom) ->
                        nom.startsWith("spins_T") && nom.endsWith(".avro"));
                if (l != null && l.length > 0) {
                    res.add(d);
                }
            }
        }
        if (res.isEmpty()) {
            // la racine elle-meme est peut-etre un dossier de run
            final String[] l = root.list((dir, nom) ->
                    nom.startsWith("spins_T") && nom.endsWith(".avro"));
            if (l != null && l.length > 0) {
                res.add(root);
            }
        }
        res.sort((a, b) -> a.getName().compareTo(b.getName()));
        return res;
    }

    private static Double temperatureDuNom(String nom) {
        final Matcher m = TAG_T.matcher(nom);
        return m.find() ? Double.valueOf(m.group(1)) : null;
    }

    private static void traiterDossier(File dossier, File inputsDir, int La, int Lb,
                                       int H, int K, int nMes, int rate, int therm,
                                       int overrelax, int threads, double sig0,
                                       long seed, String outDir,
                                       double tempCible) throws IOException {
        final long t0 = System.nanoTime();

        // ---- (1) groupe de temperature (le plus froid, ou le plus proche de --temp) ----
        final String[] snaps = dossier.list((dir, nom) ->
                nom.startsWith("spins_T") && nom.endsWith(".avro"));
        if (snaps == null || snaps.length == 0) {
            throw new SauteException("aucun spins_T*.avro");
        }
        final MapSuivi groupe = groupeLePlusFroid(dossier, snaps,
                Double.isNaN(tempCible) ? null : tempCible);
        final double T = groupe.t;

        // ---- (2) fichier d'input par nom de run ----
        // Nouveau nommage : <ddMMyyyy>[_<HHmmss>][_]... — l'underscore final est OPTIONNEL
        // (les anciens dossiers ddMMyyyy nomFichier sont collees sans separateur et doivent
        // rester parsables).
        String nomRun = dossier.getName().replaceFirst("^\\d{8}(?:_\\d{6})?_?", "");
        final Matcher mh = SUFFIXE_H.matcher(nomRun);
        double hVal = 0.0;
        if (mh.find()) {
            hVal = Double.parseDouble(mh.group(1));
            nomRun = nomRun.substring(0, mh.start());
        }
        if (hVal != 0.0) {
            throw new SauteException("champ h = " + hVal + " non supporte (outil limite a h = 0)");
        }
        final File input = new File(inputsDir, nomRun);
        if (!input.isFile()) {
            throw new SauteException("fichier d'input introuvable : " + input.getPath());
        }

        // ---- (3) reseau + maille (recette de production) ----
        final ReplicaFactory.Built built = ReplicaFactory.build(new Input(input.getPath()), La, Lb);
        final Lattice lattice = built.lattice();
        final UnitCell uc = built.uc();
        final int nSites = lattice.size;

        // ---- (4) MC relance depuis chaque instantane du groupe ; S_dimer accumule ----
        final DimerStructureFactor dimerSF = new DimerStructureFactor(H, K, lattice);
        dimerSF.setOutputDirPath(outDir);

        final double[] S = new double[3 * nSites];
        double sig = sig0;
        int nConfig = 0;
        boolean sweeperPret = false;
        CheckerboardMetropolis sweeper = null;

        for (final String snap : groupe.fichiers) {
            final File f = new File(dossier, snap);
            lireDans(f, nSites, S);                       // normalise aussi
            if (!sweeperPret) {
                // le sweeper reference le tableau du lattice : charger AVANT
                lattice.setSpin1Dlattice(S);
                sweeper = new CheckerboardMetropolis(lattice,
                        LatticeColoring.of(lattice), threads, seed);
                sweeperPret = true;
            } else {
                // meme reference que celle vue par le sweeper : copie en place
                final double[] dst = lattice.getSpin1Dlattice();
                System.arraycopy(S, 0, dst, 0, dst.length);
            }

            for (int p = 0; p < therm; p++) {
                sweeper.sweep(T, sig);
            }
            for (int m = 0; m < nMes; m++) {
                for (int p = 0; p < rate; p++) {
                    final CheckerboardMetropolis.SweepStats st = sweeper.sweep(T, sig);
                    sig = CheckerboardMetropolis.nextSigma(sig, st.rejectionRate());
                }
                for (int p = 0; p < overrelax; p++) {
                    sweeper.overRelaxationSweep();
                }
                // configuration courante, normalisee par site (comme sqwBlock statique)
                final double[] src = lattice.getSpin1Dlattice();
                for (int i = 0; i < nSites; i++) {
                    final double inv = 1.0 / MathOps.norm(src, 3 * i);
                    S[3 * i] = src[3 * i] * inv;
                    S[3 * i + 1] = src[3 * i + 1] * inv;
                    S[3 * i + 2] = src[3 * i + 2] * inv;
                }
                dimerSF.processConfig(uc, S);
                nConfig++;
            }
        }

        // ---- (5) ecriture du SEUL fichier de sortie ----
        final String tagT = String.format(java.util.Locale.US, "%.5f", T);
        final String nomFichier =
                "structure_dimere_" + nomRun + "_T" + tagT + "_h0.000000.avro";
        dimerSF.writeAvro(nomFichier);

        System.out.println(String.format(
                "OK %s : T = %.5f, %d instantane(s) x %d mesures = %d configs, "
                        + "sigma final = %.3f, %.1f s",
                nomRun, T, groupe.fichiers.length, nMes, nConfig, sig,
                (System.nanoTime() - t0) / 1e9));
        System.out.println("   -> " + nomFichier);
    }

    // ------------------------------------------------------------------

    /** Un groupe d'instantanes partageant la meme temperature (a TOL_T pres). */
    private static class MapSuivi {
        final double t;
        final String[] fichiers;
        MapSuivi(double t, String[] fichiers) { this.t = t; this.fichiers = fichiers; }
    }

    /**
     * Groupe d'instantanes a une temperature du dossier : le plus froid si
     * {@code cible == null}, sinon le groupe le plus proche de {@code cible}
     * (a TOL_T pres autour de la temperature retenue).
     */
    private static MapSuivi groupeLePlusFroid(File dossier, String[] snaps,
                                              Double cible) throws SauteException {
        double tRef = Double.POSITIVE_INFINITY;
        for (final String nom : snaps) {
            final Double t = temperatureDuNom(nom);
            if (t == null) {
                continue;
            }
            if (cible == null) {
                tRef = Math.min(tRef, t);
            } else if (tRef == Double.POSITIVE_INFINITY
                    || Math.abs(t - cible) < Math.abs(tRef - cible)) {
                tRef = t;
            }
        }
        if (tRef == Double.POSITIVE_INFINITY) {
            throw new SauteException("tag T imparsable dans les noms");
        }
        final List<String> retenus = new ArrayList<>();
        for (final String nom : snaps) {
            final Double t = temperatureDuNom(nom);
            if (t != null && t <= tRef + TOL_T) {
                retenus.add(nom);
            }
        }
        retenus.sort(String::compareTo);
        return new MapSuivi(tRef, retenus.toArray(new String[0]));
    }

    /** Lit un Avro d'instantane des spins (spins_T) dans S[3*site+g] et normalise par site. */
    private static void lireDans(File f, int nSites, double[] S) throws IOException {
        try (DataFileReader<GenericRecord> reader =
                new DataFileReader<>(f, new GenericDatumReader<GenericRecord>())) {
            int i = 0;
            while (reader.hasNext()) {
                if (i >= nSites) {
                    throw new IOException(f.getName() + " : plus de " + nSites + " sites");
                }
                final GenericRecord r = reader.next();
                S[3 * i]     = ((Number) r.get("S_x")).doubleValue();
                S[3 * i + 1] = ((Number) r.get("S_y")).doubleValue();
                S[3 * i + 2] = ((Number) r.get("S_z")).doubleValue();
                i++;
            }
            if (i != nSites) {
                throw new IOException(f.getName() + " : " + i + " sites lus, "
                        + nSites + " attendus");
            }
        }
        for (int i = 0; i < nSites; i++) {
            final double inv = 1.0 / MathOps.norm(S, 3 * i);
            S[3 * i] *= inv;
            S[3 * i + 1] *= inv;
            S[3 * i + 2] *= inv;
        }
    }
}
