import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Set;

/**
 * Coloration propre du <b>graphe d'interaction</b> d'un {@link Lattice} : deux sites i et j sont
 * adjacents si j appartient a la liste de voisins de i <em>ou</em> i a celle de j (symetrisation ;
 * les boucles i==j sont ignorees). Une coloration propre partitionne le reseau en classes de sites
 * <em>mutuellement non interagissants</em> : c'est la generalisation du damier (checkerboard) a un
 * graphe d'interaction quelconque, ce qui autorise la mise a jour Metropolis simultanee de tous les
 * sites d'une meme classe (cf. {@link CheckerboardMetropolis}).
 *
 * <h2>Etape 1 : coloration structurelle (invariante par translation)</h2>
 * Le reseau est un cristal : l'indice global se decode en (i, j, s) via
 * {@code idx = (j*La + i)*Lsr + s}. Les liens ne dependent que de l'offset
 * {@code (di, dj, s -> s2)}, dont l'ensemble est petit. On cherche alors une coloration de la forme
 * <pre>  color(i, j, s) = (alpha*i + beta*j + gamma*s) mod m </pre>
 * valide si et seulement si
 * <ol>
 *   <li>{@code (alpha*La) % m == 0} et {@code (beta*Lb) % m == 0} (coherence avec les conditions
 *       periodiques : un lien qui s'enroule doit donner la meme difference de couleur) ;</li>
 *   <li>pour tout offset, {@code (alpha*di + beta*dj + gamma*(s2 - s)) mod m != 0} (aucun lien
 *       monochrome).</li>
 * </ol>
 * On balaie m = 1..16 et alpha, beta, gamma dans 0..m-1 ; on retient le plus petit m valide
 * (departages : moins de classes non vides, puis plus petite classe maximale). Cette famille donne
 * des classes <b>parfaitement equilibrees</b> par construction des que m divise La et Lb.
 *
 * <h2>Etape 2 : repli DSATUR + reequilibrage</h2>
 * Si aucune coloration lineaire n'existe (par exemple La premier avec l'ordre de l'interaction),
 * on utilise DSATUR (Brelaz, Comm. ACM 22, 251 (1979)) sur l'adjacence symetrisee, puis un
 * reequilibrage equitable : tant que {@code max - min > 1}, on parcourt une <em>copie</em> des
 * sommets de la plus grande classe et on deplace chaque sommet vers la premiere classe c (triee par
 * taille croissante) telle que {@code taille(c) < taille(max) - 1} et qu'aucun voisin de ce sommet
 * ne porte deja la couleur c. Le test de voisinage est en O(deg) via {@code colorOf[]}.
 *
 * <p>La coloration est systematiquement validee par {@link #verify()}.</p>
 */
public final class LatticeColoring {

    /** Borne superieure du balayage structurel sur m. */
    private static final int MAX_M = 16;

    private final Lattice lattice;
    private final int[] colorOf;
    private final int[][] classes;
    private final boolean structural;
    private final String formula;

    private LatticeColoring(Lattice lattice, int[] rawColors, boolean structural, String formula) {
        this.lattice = lattice;
        this.structural = structural;
        this.formula = formula;

        // Compactage : on supprime les classes vides et on renumerote 0..C-1.
        int maxRaw = 0;
        for (int c : rawColors) {
            if (c < 0) throw new IllegalStateException("site non colorie");
            if (c > maxRaw) maxRaw = c;
        }
        int[] count = new int[maxRaw + 1];
        for (int c : rawColors) count[c]++;
        int[] remap = new int[maxRaw + 1];
        Arrays.fill(remap, -1);
        int nc = 0;
        for (int c = 0; c <= maxRaw; c++) if (count[c] > 0) remap[c] = nc++;

        this.colorOf = new int[rawColors.length];
        this.classes = new int[nc][];
        for (int c = 0; c <= maxRaw; c++) if (remap[c] >= 0) classes[remap[c]] = new int[count[c]];
        int[] fill = new int[nc];
        // Les sites sont parcourus dans l'ordre croissant : chaque classe est donc triee
        // par construction (localite spatiale des acces memoire).
        for (int site = 0; site < rawColors.length; site++) {
            int c = remap[rawColors[site]];
            colorOf[site] = c;
            classes[c][fill[c]++] = site;
        }
    }

    // ------------------------------------------------------------------ fabrique

    /**
     * Construit une coloration propre du graphe d'interaction de {@code lattice}.
     *
     * <p>La coloration structurelle est essayee en premier ; si elle existe mais echoue
     * inopinement a {@link #verify()} (offset d'interaction non couvert par la famille lineaire,
     * liste de voisins incoherente...), on <b>se replie sur DSATUR + reequilibrage</b> plutot que
     * d'echouer : le repli est toujours applicable puisqu'il travaille directement sur l'adjacence
     * symetrisee reelle. L'exception n'est relancee que si le repli echoue lui aussi ; l'echec
     * structurel y est alors joint en {@code suppressed}.</p>
     */
    public static LatticeColoring of(Lattice lattice) {
        List<int[]> offsets = collectOffsets(lattice);
        int[] best = searchStructural(lattice, offsets);
        IllegalStateException structuralFailure = null;
        if (best != null) {
            int m = best[0], alpha = best[1], beta = best[2], gamma = best[3];
            int[] raw = new int[lattice.size];
            int la = lattice.lengthA, lb = lattice.lengthB, lsr = lattice.lengthSr;
            for (int j = 0; j < lb; j++) {
                for (int i = 0; i < la; i++) {
                    for (int s = 0; s < lsr; s++) {
                        raw[(j * la + i) * lsr + s] = Math.floorMod(alpha * i + beta * j + gamma * s, m);
                    }
                }
            }
            LatticeColoring lc = new LatticeColoring(lattice, raw, true,
                    "(" + alpha + "*i + " + beta + "*j + " + gamma + "*s) mod " + m);
            try {
                lc.verify();
                return lc;
            } catch (IllegalStateException e) {
                structuralFailure = e; // repli sur DSATUR ci-dessous
            }
        }
        int[][] csr = buildAdjacency(lattice);
        int[] raw = dsatur(csr[0], csr[1]);
        rebalance(csr[0], csr[1], raw);
        LatticeColoring lc = new LatticeColoring(lattice, raw, false, "DSATUR+rebalance");
        try {
            lc.verify();
        } catch (IllegalStateException e) {
            if (structuralFailure != null) e.addSuppressed(structuralFailure);
            throw e;
        }
        return lc;
    }

    // ------------------------------------------------------------ etape 1 : structurel

    /**
     * Collecte l'ensemble deduplique des offsets d'interaction {@code (di, dj, s, s2)} sur
     * <em>tous</em> les sites. Les deux sens sont ajoutes (symetrisation explicite) ; les boucles
     * sont ignorees.
     */
    private static List<int[]> collectOffsets(Lattice lattice) {
        final int la = lattice.lengthA, lb = lattice.lengthB, lsr = lattice.lengthSr;
        final List<int[]> out = new ArrayList<>();
        final Set<Long> seen = new HashSet<>();
        for (int site = 0; site < lattice.size; site++) {
            int s = site % lsr;
            int rest = site / lsr;
            int i = rest % la;
            int j = rest / la;
            int[] nb = lattice.getNeighborSites(site);
            if (nb == null) continue;
            for (int k = 0; k < nb.length; k++) {
                int other = nb[k];
                if (other == site) continue; // boucle : ignoree
                int s2 = other % lsr;
                int rest2 = other / lsr;
                int i2 = rest2 % la;
                int j2 = rest2 / la;
                int di = wrap(i2 - i, la);
                int dj = wrap(j2 - j, lb);
                addOffset(out, seen, di, dj, s, s2);
                addOffset(out, seen, -di, -dj, s2, s); // arete inverse (symetrisation)
            }
        }
        return out;
    }

    private static void addOffset(List<int[]> out, Set<Long> seen, int di, int dj, int s, int s2) {
        if (di == 0 && dj == 0 && s == s2) return;
        long key = ((((long) (di + (1 << 20)) << 21) | (dj + (1 << 20))) << 16) | ((long) s << 8) | s2;
        if (seen.add(key)) out.add(new int[] { di, dj, s, s2 });
    }

    /** Replie d dans (-L/2, L/2]. */
    private static int wrap(int d, int len) {
        int r = Math.floorMod(d, len);
        if (r > len / 2) r -= len;
        return r;
    }

    /**
     * Balaye m = 1..16 puis alpha, beta, gamma dans 0..m-1.
     *
     * @return {m, alpha, beta, gamma} du meilleur candidat, ou null si la famille lineaire echoue.
     */
    private static int[] searchStructural(Lattice lattice, List<int[]> offsets) {
        final int la = lattice.lengthA, lb = lattice.lengthB, lsr = lattice.lengthSr;
        for (int m = 1; m <= MAX_M; m++) {
            int[] best = null;
            int bestNonEmpty = Integer.MAX_VALUE;
            long bestMax = Long.MAX_VALUE;
            for (int alpha = 0; alpha < m; alpha++) {
                if ((alpha * la) % m != 0) continue;              // coherence periodique en a
                for (int beta = 0; beta < m; beta++) {
                    if ((beta * lb) % m != 0) continue;           // coherence periodique en b
                    for (int gamma = 0; gamma < m; gamma++) {
                        if (!isValid(offsets, m, alpha, beta, gamma)) continue;
                        long[] sizes = classSizes(m, alpha, beta, gamma, la, lb, lsr);
                        int nonEmpty = 0;
                        long maxSize = 0;
                        for (long v : sizes) {
                            if (v > 0) nonEmpty++;
                            if (v > maxSize) maxSize = v;
                        }
                        if (nonEmpty < bestNonEmpty || (nonEmpty == bestNonEmpty && maxSize < bestMax)) {
                            bestNonEmpty = nonEmpty;
                            bestMax = maxSize;
                            best = new int[] { m, alpha, beta, gamma };
                        }
                    }
                }
            }
            if (best != null) return best; // plus petit m valide
        }
        return null;
    }

    private static boolean isValid(List<int[]> offsets, int m, int alpha, int beta, int gamma) {
        for (int k = 0; k < offsets.size(); k++) {
            int[] o = offsets.get(k);
            int d = alpha * o[0] + beta * o[1] + gamma * (o[3] - o[2]);
            if (Math.floorMod(d, m) == 0) return false;
        }
        return true;
    }

    /** Tailles des m classes, calculees par convolution des histogrammes (independant de N). */
    private static long[] classSizes(int m, int alpha, int beta, int gamma, int la, int lb, int lsr) {
        long[] hi = new long[m], hj = new long[m], hs = new long[m];
        for (int i = 0; i < la; i++) hi[Math.floorMod(alpha * i, m)]++;
        for (int j = 0; j < lb; j++) hj[Math.floorMod(beta * j, m)]++;
        for (int s = 0; s < lsr; s++) hs[Math.floorMod(gamma * s, m)]++;
        long[] ij = new long[m];
        for (int a = 0; a < m; a++) {
            if (hi[a] == 0) continue;
            for (int b = 0; b < m; b++) ij[(a + b) % m] += hi[a] * hj[b];
        }
        long[] out = new long[m];
        for (int a = 0; a < m; a++) {
            if (ij[a] == 0) continue;
            for (int b = 0; b < m; b++) out[(a + b) % m] += ij[a] * hs[b];
        }
        return out;
    }

    // ------------------------------------------------------------ etape 2 : DSATUR

    /** Adjacence symetrisee, dedupliquee, au format CSR : {rowStart[n+1], adj[]}. */
    private static int[][] buildAdjacency(Lattice lattice) {
        final int n = lattice.size;
        int[] deg = new int[n];
        for (int i = 0; i < n; i++) {
            int[] nb = lattice.getNeighborSites(i);
            if (nb == null) continue;
            for (int k = 0; k < nb.length; k++) {
                int j = nb[k];
                if (j == i) continue;
                deg[i]++;
                deg[j]++;
            }
        }
        int[] start = new int[n + 1];
        for (int i = 0; i < n; i++) start[i + 1] = start[i] + deg[i];
        int[] adj = new int[start[n]];
        int[] pos = Arrays.copyOf(start, n);
        for (int i = 0; i < n; i++) {
            int[] nb = lattice.getNeighborSites(i);
            if (nb == null) continue;
            for (int k = 0; k < nb.length; k++) {
                int j = nb[k];
                if (j == i) continue;
                adj[pos[i]++] = j;
                adj[pos[j]++] = i;
            }
        }
        // deduplication ligne par ligne
        int[] outStart = new int[n + 1];
        int write = 0;
        for (int i = 0; i < n; i++) {
            int from = start[i], to = start[i + 1];
            Arrays.sort(adj, from, to);
            outStart[i] = write;
            int prev = -1;
            for (int k = from; k < to; k++) {
                if (adj[k] != prev) {
                    adj[write++] = adj[k];
                    prev = adj[k];
                }
            }
        }
        outStart[n] = write;
        return new int[][] { outStart, Arrays.copyOf(adj, write) };
    }

    /** DSATUR (Brelaz 1979) : on colorie d'abord le sommet de saturation maximale. */
    private static int[] dsatur(int[] rowStart, int[] adj) {
        final int n = rowStart.length - 1;
        final int[] color = new int[n];
        Arrays.fill(color, -1);
        int maxDeg = 0;
        for (int i = 0; i < n; i++) maxDeg = Math.max(maxDeg, rowStart[i + 1] - rowStart[i]);
        final int words = Math.max(1, (maxDeg + 2 + 63) >>> 6);
        final long[] used = new long[n * words];
        final int[] sat = new int[n];

        // File de priorite avec suppression paresseuse : (saturation desc, degre desc, indice asc).
        PriorityQueue<int[]> pq = new PriorityQueue<>((x, y) -> {
            if (x[0] != y[0]) return y[0] - x[0];
            if (x[1] != y[1]) return y[1] - x[1];
            return x[2] - y[2];
        });
        for (int i = 0; i < n; i++) pq.add(new int[] { 0, rowStart[i + 1] - rowStart[i], i });

        int colored = 0;
        while (colored < n) {
            int[] e = pq.poll();
            if (e == null) break;
            int v = e[2];
            if (color[v] >= 0 || sat[v] != e[0]) continue; // entree perimee
            int base = v * words;
            int c = 0;
            while (c <= maxDeg + 1 && (used[base + (c >>> 6)] & (1L << (c & 63))) != 0L) c++;
            color[v] = c;
            colored++;
            for (int k = rowStart[v]; k < rowStart[v + 1]; k++) {
                int u = adj[k];
                if (color[u] >= 0) continue;
                int ub = u * words + (c >>> 6);
                long bit = 1L << (c & 63);
                if ((used[ub] & bit) == 0L) {
                    used[ub] |= bit;
                    sat[u]++;
                    pq.add(new int[] { sat[u], rowStart[u + 1] - rowStart[u], u });
                }
            }
        }
        return color;
    }

    /**
     * Reequilibrage equitable des classes. Corrige les deux defauts de la version naive :
     * on itere sur une <em>copie</em> des membres de la classe (jamais sur la liste mutee) et on
     * essaie <em>toutes</em> les classes triees par taille croissante, pas seulement la plus petite.
     */
    private static void rebalance(int[] rowStart, int[] adj, int[] color) {
        int nc = 0;
        for (int c : color) nc = Math.max(nc, c + 1);
        if (nc <= 1) return;
        final int n = color.length;
        int[] size = new int[nc];
        for (int c : color) size[c]++;

        final int maxPasses = 64;
        for (int pass = 0; pass < maxPasses; pass++) {
            int largest = 0, smallest = 0;
            for (int c = 1; c < nc; c++) {
                if (size[c] > size[largest]) largest = c;
                if (size[c] < size[smallest]) smallest = c;
            }
            if (size[largest] - size[smallest] <= 1) return;

            // Copie instantanee des membres de la plus grande classe.
            int[] members = new int[size[largest]];
            int p = 0;
            for (int v = 0; v < n; v++) if (color[v] == largest) members[p++] = v;

            // Classes candidates triees par taille croissante.
            Integer[] order = new Integer[nc];
            for (int c = 0; c < nc; c++) order[c] = c;
            final int[] sizeSnapshot = size.clone();
            Arrays.sort(order, (x, y) -> Integer.compare(sizeSnapshot[x], sizeSnapshot[y]));

            boolean moved = false;
            for (int idx = 0; idx < members.length; idx++) {
                int v = members[idx];
                if (color[v] != largest) continue;
                for (int oi = 0; oi < nc; oi++) {
                    int c = order[oi];
                    if (c == largest) continue;
                    if (size[c] >= size[largest] - 1) continue; // taille reevaluee dynamiquement
                    if (hasNeighborColored(rowStart, adj, color, v, c)) continue;
                    color[v] = c;
                    size[c]++;
                    size[largest]--;
                    moved = true;
                    break;
                }
            }
            if (!moved) return; // passe complete sans mouvement
        }
    }

    private static boolean hasNeighborColored(int[] rowStart, int[] adj, int[] color, int v, int c) {
        for (int k = rowStart[v]; k < rowStart[v + 1]; k++) if (color[adj[k]] == c) return true;
        return false;
    }

    // ------------------------------------------------------------------- API

    /**
     * Reseau colorie. Permet de verifier qu'une coloration est bien utilisee avec le reseau dont
     * elle provient (cf. le constructeur de {@link CheckerboardMetropolis}).
     */
    public Lattice lattice() { return lattice; }

    /** Nombre de classes de couleur non vides. */
    public int numColors() { return classes.length; }

    /** Couleur du site donne, dans 0..numColors()-1. */
    public int colorOf(int site) { return colorOf[site]; }

    /**
     * Sites de la classe {@code c}, tries par ordre croissant.
     * <b>Le tableau retourne est interne : il ne doit pas etre modifie.</b>
     */
    public int[] classSites(int c) { return classes[c]; }

    /** Toutes les classes. <b>Tableau interne : ne pas modifier.</b> */
    public int[][] classes() { return classes; }

    /** true si la coloration provient de la famille lineaire (equilibree par construction). */
    public boolean isStructural() { return structural; }

    /** Formule de coloration, par exemple "(1*i + 1*j + 2*s) mod 6" ou "DSATUR+rebalance". */
    public String formula() { return formula; }

    /** Taille de la plus grande classe. */
    public int maxClassSize() {
        int m = 0;
        for (int[] cl : classes) m = Math.max(m, cl.length);
        return m;
    }

    /** Taille de la plus petite classe. */
    public int minClassSize() {
        int m = Integer.MAX_VALUE;
        for (int[] cl : classes) m = Math.min(m, cl.length);
        return classes.length == 0 ? 0 : m;
    }

    /** Resume sur une ligne. */
    public String summary() {
        StringBuilder sb = new StringBuilder();
        sb.append("LatticeColoring[").append(numColors()).append(" colors, ")
          .append(structural ? "structural " : "").append(formula).append(", sizes=");
        for (int c = 0; c < classes.length; c++) {
            if (c > 0) sb.append('/');
            sb.append(classes[c].length);
        }
        sb.append(", balanced=").append(maxClassSize() - minClassSize() <= 1).append(']');
        return sb.toString();
    }

    /**
     * Verifie qu'aucune arete n'est monochrome et que tous les sites sont colories.
     *
     * @throws IllegalStateException si la coloration n'est pas propre
     */
    public void verify() {
        for (int site = 0; site < lattice.size; site++) {
            int c = colorOf[site];
            if (c < 0 || c >= classes.length) {
                throw new IllegalStateException("site " + site + " non colorie (couleur " + c + ")");
            }
            int[] nb = lattice.getNeighborSites(site);
            if (nb == null) continue;
            for (int k = 0; k < nb.length; k++) {
                int other = nb[k];
                if (other == site) continue;
                if (colorOf[other] == c) {
                    throw new IllegalStateException(
                            "arete monochrome " + site + "-" + other + " (couleur " + c + ")");
                }
            }
        }
    }

    /**
     * Test de proprete autonome, utilisable sur un tableau de couleurs quelconque
     * (pratique pour les tests : on peut lui passer une coloration volontairement fausse).
     */
    public static boolean isProper(Lattice lattice, int[] colorOf) {
        if (colorOf == null || colorOf.length < lattice.size) return false;
        for (int site = 0; site < lattice.size; site++) {
            if (colorOf[site] < 0) return false;
            int[] nb = lattice.getNeighborSites(site);
            if (nb == null) continue;
            for (int k = 0; k < nb.length; k++) {
                int other = nb[k];
                if (other == site) continue;
                if (colorOf[other] == colorOf[site]) return false;
            }
        }
        return true;
    }
}
