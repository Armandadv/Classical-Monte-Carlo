# Monte Carlo Classique

Code de simulation Monte Carlo classique de spins que j'ai développé durant ma thèse puis poursuivit après. Il
calcule l'état fondamental d'un Hamiltonien de spins anisotropes sur réseau en nid d'abeille,
les propriétés thermodynamiques ainsi que les excitations associées. La méthode permet de
calculer les observables physiques, telles que l'aimantation, le
facteur de structure magnétique statique et dynamique, ou encore la susceptibilité en champ
alternatif. Ces observables sont ensuite comparées aux données expérimentales pour choisir
parmi les modèles ceux qui sont les plus prometteurs.

Le code est écrit en Java (JDK 21) et ne dépend que des bibliothèques fournies dans `lib/`
(Avro pour les sorties, JTransforms pour les FFT, Commons Math pour les opérations arithémtiques).


## Le modèle simulé

Le code simule un Hamiltonien classique général d'interactions anisotropes. Une solution
consiste à généraliser le terme d'échange de Heisenberg, qui dans sa forme la plus simple
est un scalaire, par une matrice de dimension 3 qui sera notée $\mathbb{J}$. De cette
manière, le choix des éléments de cette matrice permet de favoriser une certaine direction
et ainsi d'introduire de l'anisotropie d'échange. L'anisotropie de l'ion libre est portée par
un tenseur $g$ : appliquée à des spins aux sites $i$ et $j$, elle s'écrit comme
$\mathbf{S}_i = g\,\mathbf{s}_i$ et $\mathbf{S}_j = g\,\mathbf{s}_j$, soit
$\mathbf{S}_i^{\,T}\mathbb{J}\mathbf{S}_j = \mathbf{s}_i^{\,T}\left(g^T\mathbb{J}g\right)\mathbf{s}_j$,
c'est-à-dire qu'elle est directement incorporée dans les termes d'échange anisotropes. Le
terme de Zeeman en revanche couple le champ appliqué au spin habillé par $g$. Le modèle
simulé s'écrit donc :

$$H = -\sum_{n}\ \sum_{\langle ij\rangle_n} \mathbf{S}_i \cdot \mathbb{J}_{n} \cdot \mathbf{S}_j
\;-\; \sum_i \mathbf{B}\cdot g\cdot\mathbf{S}_i
\;-\; b\sum_{\langle ij\rangle\in\gamma} \left(S_i^\gamma S_j^\gamma\right)^2$$

Le premier terme porte sur les $n$èmes voisins, chaque $\mathbb{J}_n$ étant une matrice
symétrique de dimension 3 lue dans le fichier d'entrée, la symétrie étant contrôlée à la
construction et la réciprocité des liens vérifiée et signalée. Toute paramétrisation usuelle
s'y exprime pour des premiers voisins edans la base du fichier d'input. Par exemple dans la base Kitaev, on peut écrire la matrice qui couple le liens $z$ comme :

$$\mathbb{J}_z = \begin{pmatrix} J_1 & \Gamma & \Gamma'\\ \Gamma & J_1 & \Gamma'\\ \Gamma' & \Gamma' & J_1 + K \end{pmatrix},$$

où le terme commun sur la diagonale $J_1$ correspond au couplage isotrope de Heisenberg, le
terme $K$ est le terme de Kitaev qui couple les composantes de spin le long du lien
considéré, et $\Gamma$, $\Gamma'$ sont les échanges hors diagonaux. Les voisins plus lointains
sont peuvent eux aussi être écrit sous cette forme $\mathbb{J}_2$, $\mathbb{J}_3$, $\mathbb{J}_4$. Le dernier terme du Hamiltonien, optionnel
(activé par une ligne `QUARTIC` dans le fichier d'entrée), est un couplage biquadratique sur
premiers voisins le long du lien $\gamma$. Sa variable de liaison $\tau_b =
(S_i^\gamma S_j^\gamma)^2$ est mesurée par un canal de sortie dédié.

Les conventions sont les suivantes : 

- $k_B = 1$ : la température est dans les mêmes unités que les couplages $J$ du fichier d'entrée (meV).
- Spins unitaires : les spins sont des vecteurs unitaires 3D ($|\mathbf{S}| = 1$).
- $\hbar = 1$ dans la dynamique.
- Champ magnétique : le champ $\mathbf{B}$ (en meV) est exprimé via le tenseur $g$ du fichier d'entrée.
- Conditions aux bords périodiques : elles imposent une quantification des vecteurs d'onde ($\mathbf{Q} = 2\pi\mathbf{n}/L$) et provoquent un verrouillage commensurable intrinsèque.
- Tailles de réseau : typiquement de $24\times24\times2$ pour les mesures sous champ magnétique, et jusqu'à $60\times60\times2$ lorsque l'état d'équilibre est supposé incommensurable.

## Échantillonnage markovien et algorithme de Metropolis

Les méthodes de Monte Carlo sont une famille de simulations stochastiques permettant
d'estimer numériquement des intégrales. Elles sont particulièrement efficaces pour
échantillonner l'espace des phases à partir d'une distribution de probabilité, même lorsque
celle-ci n'est connue qu'à une constante de normalisation près. Dans le cas de la physique
statistique, elles sont souvent utilisées pour calculer les moyennes à l'équilibre
thermodynamique d'une observable physique telle que l'aimantation, l'énergie ou l'entropie
en se basant sur la distribution de Boltzmann, où la probabilité d'occuper un micro-état
$x$ est donnée par :

$$p(x) = \frac{1}{Z}\, e^{-\beta E(x)}.$$

L'idée de Metropolis est de construire une chaîne de Markov dont la distribution stationnaire
est précisément $p(x)$. Une condition suffisante pour que $p(x)$ soit stationnaire est de
respecter la condition de bilan détaillé (*detailed balance*) :

$$p(x)\,K_{x,x'} = p(x')\,K_{x',x} \qquad \forall\, x, x' \in \Omega.$$

Dans le cas d'une distribution de proposition symétrique et de l'utilisation du critère de
Metropolis, la condition d'acceptation s'écrit :

$$\alpha(x,x') = \min\left\{1,\ e^{-\beta\left(E(x')-E(x)\right)}\right\}.$$

La condition d'acceptation ou de rejet est déterminée en tirant un nombre aléatoire
$r\in[0,1]$ dans une distribution uniforme, et en comparant $r$ au terme
$e^{-\beta(E(x')-E(x))}$. Cela permet au système de rester en mouvement, de ne pas rester
dans un minimum local et d'explorer l'espace des configurations qui semblent moins
favorables. Ainsi, la probabilité d'acceptation ne dépend pas de la fonction de partition.
Un pas Monte Carlo (MCS) correspond à $N$ tentatives de mise à jour des spins par
l'algorithme de Metropolis sur un système de $N$ spins, de sorte qu'en moyenne pour un pas
chaque spin est testé une fois. Le générateur de nombres pseudo-aléatoires est un
Xoshiro256++.

### Mouvement gaussien

Dans le cas des systèmes possédant des spins continus, l'échantillonnage de l'espace des
phases ne peut pas se faire simplement en retournant le spin comme on le ferait dans un
modèle d'Ising ou de Potts : les spins peuvent varier de manière infinitésimale et prendre
n'importe quelle valeur sur une sphère. Une manière triviale pour l'échantillonner serait de
choisir de manière aléatoire une nouvelle orientation de spin sur une sphère. Ce choix est
efficace dans le régime paramagnétique, mais moins pertinent à basse température, où l'on
s'attend plutôt à ce que les spins fluctuent autour de leur champ moléculaire via de faibles
déviations. Pour surmonter cette difficulté, la méthode du mouvement gaussien consiste à
adapter pour chaque température l'angle solide d'échantillonnage pour maintenir un taux
d'acceptation cible. La nouvelle position de spin est générée en prenant la somme du spin
actuel et d'un vecteur aléatoire $\Gamma$ suivant une loi normale, pondérée par un paramètre
$\sigma$ :

$$\mathbf{S}_{new} = \frac{\mathbf{S}_i + \sigma\Gamma}{|\mathbf{S}_i + \sigma\Gamma|},
\qquad \sigma_{new} = \frac{0.5}{1 - R}\,\sigma,$$

avec $R$ le taux d'acceptation local calculé à chaque pas. La valeur de $\sigma$ est ajustée
dynamiquement pour maintenir un taux d'acceptation à 50 % à chaque température. Le caractère
gaussien contenu dans la variable $\sigma$ permet d'explorer au-delà d'un angle limite
imposé : cette méthode permet d'explorer plus efficacement, évite de se piéger dans certains
minima locaux et garantit un échantillonnage efficace à toutes les températures, ce qui est
parfaitement adapté pour une simulation du recuit simulé. Pendant la production en revanche,
$\sigma$ est figé à sa valeur acquise en thermalisation : l'adapter en cours de mesure
rendrait la proposition dépendante de l'historique des taux de rejet, donc la chaîne non
markovienne.

## Sur-relaxation et dynamique moléculaire de spins

Malgré ses avantages, la méthode d'échantillonnage par mouvement gaussien présente des
limitations intrinsèques à la physique des systèmes magnétiques sur des grandes tailles de
réseau. Dans la majorité des cas, le système va générer des parois de domaines qui se
traduisent par l'augmentation de configurations métastables proches en énergie. Une même
paroi pouvant se déplacer n'importe où dans le système sans coûter plus d'énergie. Afin de
réduire ce phénomène dans les simulations et d'augmenter l'ergodicité, nous avons combiné
la méthode d'échantillonnage avec de la dynamique. La première composante est une
sur-relaxation microcanonique (Creutz 1987, Brown et Woch 1987). Le champ moléculaire local
ressenti par le spin $i$ s'écrit :

$$\mathbf{h}_i = \sum_{k\in\text{voisins}} \mathbb{J}_{ik}\,\mathbf{S}_k + g^T\mathbf{B},$$

et la sur-relaxation consiste à remplacer le spin par son reflet par rapport à la direction
de ce champ :

$$\mathbf{S}_{new} = 2\left(\mathbf{S}_{old}\cdot\hat{h}_i\right)\hat{h}_i - \mathbf{S}_{old}.$$

La projection du spin sur son champ est inchangée par cette réflexion,
$\mathbf{S}_{new}\cdot\mathbf{h}_i = \mathbf{S}_{old}\cdot\mathbf{h}_i$, de sorte que l'énergie
locale $-\mathbf{S}\cdot\mathbf{h}_i$ est exactement conservée, les voisins étant figés
pendant l'opération. Le mouvement est donc toujours accepté : il décorrèle le système sans
aucun tirage aléatoire ni rejet, et coûte bien moins qu'une tentative Metropolis. Elle est
entrelacée avec les balayages Metropolis (deux passes par balayage par défaut) et interdite
si le terme quartique est actif, la réflexion ne conservant que le terme quadratique de
l'énergie.

La deuxième méthode s'appuie sur la résolution numérique des équations du mouvement
décrites par les équations de Landau-Lifshitz qui régissent l'évolution temporelle de chaque
spin autour de son champ moléculaire :

$$\frac{d\mathbf{S}_i(t)}{dt} = \mathbf{S}_i(t) \times \frac{\partial H}{\partial \mathbf{S}_i(t)},
\qquad \forall i \in \{1,\ldots,N\}.$$

Cette méthode est complémentaire à l'échantillonnage gaussien et permet de décorréler le
système d'une autre manière que par l'algorithme de Metropolis. De plus, comme ces équations
conservent l'énergie globale du système (pas de terme d'amortissement), cela n'implique
aucun rejet, tout en améliorant la mobilité des spins. Cette relaxation est appliquée en fin
de chaque température. Outre l'amélioration de l'ergodicité, elle permet d'accéder à une
échelle de temps réel, ce que ne permet pas la méthode Metropolis : c'est elle qui est
utilisée pour calculer le facteur de structure dynamique.

## Recuit simulé et thermalisation automatique

Pour obtenir un état d'équilibre proche de la distribution stationnaire, les simulations
sont initialisées dans un régime de haute température (régime paramagnétique) pour ensuite
être refroidies jusqu'à la température souhaitée : c'est le principe du recuit simulé. Comme
le critère d'acceptation suit une loi de Boltzmann, lorsque la température diminue, $\beta$
augmente et le taux d'acceptation tend vers zéro. Les nouvelles configurations ne sont alors
que très rarement acceptées et le système peut se figer dans une configuration métastable.
Pour laisser le temps au système de se thermaliser, on impose aux valeurs de la température
de suivre une suite géométrique $T_{i+1} = k\,T_i$ avec $k = 0.95$. Un refroidissement
adaptatif est également disponible, selon le critère de Huang, Romeo et
Sangiovanni-Vincentelli (1986) : l'amplitude des fluctuations d'énergie
$\sigma_E = \sqrt{\mathrm{Var}(E)} = T\sqrt{Nc_v}$, où $c_v = \mathrm{Var}(E)/(NT^2)$ est
la chaleur spécifique par site, prescrit le pas en température,

$$\frac{T_{k+1}}{T_k} = \exp\left(-\frac{\lambda}{\sqrt{c_v}}\right),$$

avec $\lambda$ un paramètre adimensionné : les pas sont resserrés là où la chaleur spécifique
est grande, au voisinage des transitions, et relâchés là où les fluctuations d'énergie sont
faibles.

À chaque pas en température, il est nécessaire de réaliser un certain nombre de pas Monte
Carlo qui serviront à thermaliser le système avant de mesurer les valeurs moyennes. Plutôt
qu'un nombre fixé à l'avance, le code détecte l'équilibre. Le critère principal est un test
de Geweke appliqué à deux observables, l'énergie et la norme de l'aimantation. La série
mesurée depuis le début de la température est coupée en deux fenêtres : une fenêtre
initiale $A$, qui contient la fraction la plus ancienne de la chaîne, et une fenêtre finale
$B$, qui contient la fraction la plus récente. Le test compare les moyennes de l'observable
sur les deux fenêtres par le score

$$z = \frac{\langle\mathcal{O}\rangle_A - \langle\mathcal{O}\rangle_B}
{\sqrt{v_A\, g_A/n_A + v_B\, g_B/n_B}},$$

où $\mathcal{O}$ désigne l'observable suivie (l'énergie ou la norme de l'aimantation),
$\langle\mathcal{O}\rangle_A$ et $\langle\mathcal{O}\rangle_B$ ses moyennes sur les fenêtres
initiale et finale, $v_A$ et $v_B$ ses variances dans chacune des deux fenêtres, $g_A$ et
$g_B$ les inefficacités statistiques des fenêtres, et $n_A$ et $n_B$ leurs nombres
d'échantillons.  Sous
l'hypothèse de stationnarité, $z$ suit approximativement une loi normale centrée réduite,
la thermalisation s'étend jusqu'à ce que $|z|$ repasse sous le seuil pour les deux
observables.

La production qui suit vise un nombre effectif d'échantillons indépendants plutôt qu'un
nombre de pas fixe. Si $\tau$ est le temps d'autocorrélation de l'énergie, compté en pas
Monte Carlo, l'inefficacité statistique vaut $g = 1 + 2\tau$ : c'est le facteur par lequel
il faut diviser le nombre de mesures $N_t$ pour obtenir le nombre de mesures réellement
indépendantes,

$$N_{eff} = \frac{N_t}{g}.$$

Les mesures sont espacées d'environ deux temps de corrélation, et la production est étendue
automatiquement jusqu'à livrer le $N_{eff}$ promis, le manque étant converti en balayages
supplémentaires par le facteur $g$. Un détecteur de stagnation peut en outre, si l'arrêt
anticipé correspondant est activé, mettre fin à la thermalisation ou à la production
lorsqu'il constate que le gain de $N_{eff}$ sature.

## Observables et facteur de structure dynamique

Après avoir atteint l'équilibre, on commence la mesure des quantités thermodynamiques.
L'aimantation par spin suivant un vecteur $w$ est calculée par :

$$\langle m_w\rangle = \frac{1}{N N_t}\sum_i^N \sum_n^{N_t} \mathbf{S}_i(\tau_n)\cdot w,$$

et de même pour l'énergie :

$$\langle E\rangle = -\frac{1}{2NN_t}\sum_i^N \sum_n^{N_t}
\mathbf{S}_i(\tau_n) \sum_{k\in\text{voisins de } i} \mathbb{J}_{ik}\,\mathbf{S}_k(\tau_n),$$

avec la somme sur $n$ jusqu'à $N_t$ qui représente le nombre de fois où l'on répète la
mesure pour une température donnée. Toutes les mesures, leurs erreurs et les diagnostics de
corrélation sont écrits, une ligne par température, dans `equilibration.csv`.

Une fois à la température souhaitée et à l'équilibre thermodynamique, on intègre la
dynamique des spins à partir de la résolution numérique des équations de Landau-Lifshitz
pour obtenir l'évolution temporelle de chaque spin. On calcule ensuite la double transformée
de Fourier dans le temps et dans l'espace de chaque spin :

$$\tilde{S}_\perp(Q,\omega) = \frac{1}{NN_t}\sum_{\alpha\in\{1,2\}}\sum_j^{N/2}\sum_n^{N_t}
\mathbf{S}_{j,\alpha}^{\perp}(t_n)\, e^{-iQ\cdot(R_j + r_\alpha)}\, e^{-i\omega t_n},$$

où $\mathbf{S}_{j,\alpha}^{\perp}(t_n)$ est la projection du spin habillé par le tenseur
$g$, perpendiculaire au vecteur de diffusion $Q$, les deux sous-réseaux du nid d'abeille
étant combinés par le facteur de phase $e^{-iQ\cdot(r_B-r_A)}$. Le facteur de structure dynamique est calculé via
l'expression :

$$\langle S(Q,\omega)\rangle = \tilde{S}_\perp(Q,\omega)\cdot\tilde{S}_\perp^*(Q,\omega),$$

où $\langle\cdot\rangle$ désigne la moyenne effectuée sur les instantanés successifs, le
spectre étant accumulé en moyenne courante puis écrit au dernier instantané. Il en découle
l'intensité en unité arbitraire :

$$I(Q,\omega) \propto |f(Q)|^2\, S(Q,\omega),$$

avec $f(Q)$ le facteur de forme magnétique de l'ion. Le nombre d'instantanés, la durée
d'intégration et la résolution en énergie sont réglables. Une fenêtre gaussienne
d'apodisation est appliquée aux séries temporelles avant la transformée de Fourier. Pour un vecteur de
diffusion strictement bidimensionnel ($q_z = 0$), l'hermiticité de la transformée permet de
ne calculer que la moitié des blocs de Fourier et de remplir le reste par symétrie miroir.
La composante $q_z$ est par ailleurs prise en compte dans le facteur de phase, et
enregistrée
dans le schéma des fichiers de sortie.

Deux canaux de structure supplémentaires accompagnent le canal dipolaire : le canal couleur,
qui mesure $n_i^a = (S_i^a)^2$ dans la base de Kitaev, et le canal dimère, qui mesure la
variable de liaison $\tau_b = (S_i^\gamma S_j^\gamma)^2$ du terme quartique aux milieux de
liaisons. Ces deux canaux ne sont calculés et écrits que si le terme quartique est actif.

## Susceptibilité en champ alternatif

Le code calcule également la susceptibilité AC. Le champ magnétique imposé au système
dépend du temps :

$$\mathbf{H}(t) = \mathbf{H}_{stat} + H_0 \sin(\omega t)\; g\,\hat{u},$$

où $H_0$ est l'amplitude du champ oscillant (l'option `--acH0`) en meV avant
habillage par le tenseur $g$, $\mathbf{H}_{stat}$ est le champ statique de la réplique
construit à partir du champ du balayage en champ fixe et $\hat{u}$ est la direction
d'application unitaire choisie (l'option `--acFieldDir`. La pulsation
$\omega = 2\pi/P$ est fixée par la période $P$ de l'excitation, comptée en balayages. Pour
chaque période d'excitation, l'échelle de refroidissement complète est rejouée depuis la plus haute
température de la grille : le système est d'abord équilibré sous le seul champ statique,
puis la mesure s'effectue sous le champ complet, à l'issue du régime transitoire. L'aimantation
en phase et en quadrature avec l'excitation est démodulée sur un nombre entier de périodes :

$$\chi' = \frac{1}{n\,H_0\,N}\sum_{j=1}^{n} m(t_j)\sin(\omega t_j),
\qquad \chi'' = -\frac{1}{n\,H_0\,N}\sum_{j=1}^{n} m(t_j)\cos(\omega t_j),$$

où $m(t_j) = \sum_i^N \mathbf{S}_i(t_j)\cdot\hat{u}_g$ est l'aimantation totale du système
projetée sur la direction d'excitation habillée et normalisée,
$\hat{u}_g = g\,\hat{u}/|g\,\hat{u}|$, $N$ est le nombre de
spins, la susceptibilité étant ainsi une grandeur intensive et où $n$ est le nombre de
mesures accumulées pendant la fenêtre de mesure. La somme sur un nombre entier de périodes
assure l'orthogonalité exacte des
projections sur sinus et cosinus. Le signe est choisi pour que la partie dissipative
$\chi''$ soit positive pour une réponse retardée, comme sur une mesure expérimentale.
L'horloge de la mesure est la dynamique stochastique des balayages,
chaque balayage étant une tentative de relaxation thermiquement activée, et non l'échelle de
temps réel de la précession. La sur-relaxation, active pendant l'équilibrage, est coupée
pendant la mesure. Les résultats sont moyennés sur plusieurs répliques indépendantes et
écrits dans les fichiers de sortie.

## Parallélisation : décomposition en réseau de damier

Dans l'algorithme de Metropolis, on échantillonne l'espace des phases d'un seul spin à la
fois (*single spin flip*), ce qui entraîne des temps de calcul longs sur des réseaux de spins
de grande taille. Paralléliser naïvement l'algorithme de Metropolis poserait un problème :
cela ne respecterait plus la condition de bilan détaillé. En effet, si deux spins voisins
sont sélectionnés pour être mis à jour, la modification de l'un affectera le calcul de
l'énergie de l'autre, faussant ainsi la probabilité d'acceptation $\alpha(x,x')$. Pour
contourner ce problème, la décomposition en réseau de damier (*checkerboard decomposition*)
permet de regrouper les spins qui peuvent être mis à jour indépendamment : on décompose le
réseau en sous-réseaux $\{n_0, n_1, \cdots\}$ tels que les spins d'un même sous-réseau
n'interagissent qu'avec des spins d'autres sous-réseaux. Puisque les champs locaux ne
dépendent alors que de spins figés pendant la phase, la mise à jour simultanée de tous les
spins d'un même sous-réseau est rigoureusement équivalente à un balayage séquentiel et
laisse invariante la mesure de Boltzmann.

La décomposition est d'abord cherchée sous forme structurelle : une coloration linéaire
$(\alpha i + \beta j + \gamma s) \bmod m$, au nombre minimal de couleurs, dont les classes
sont équilibrées par construction. Si aucune coloration linéaire n'existe, la décomposition
est obtenue par coloration de graphe : un algorithme de coloration (DSATUR) assigne les
couleurs en veillant à ce qu'aucun site adjacent ne partage la même couleur, puis un
rééquilibrage transfère des sites des plus grandes classes vers les plus petites jusqu'à un
écart d'au plus un site. Cette formulation permet de l'adapter à tout type de réseau (carré,
kagome, triangulaire) et à des interactions allant jusqu'aux voisins lointains. À
l'intérieur d'une réplique, les sous-réseaux sont mis à jour en parallèle par un groupe de
threads. 

## Utilisation

### Compilation

```bash
./build.sh          # compile src/ -> bin/
```
Windows : le séparateur de classpath est `;`. Linux et macOS : `:`. Sur machine partagée,
épingler la JVM avec `-XX:ActiveProcessorCount=N`.
```bash
java -cp "bin;lib/*" AppParallel ...  # sous Windows
```
```bash
java -cp "bin:lib/*" AppParallel ...     # sous Linux ou macOS
```

### Arguments positionnels

```bash
java -cp "bin:lib/*" AppParallel <inputDir> <outputDir> [totalThreads] [threadsPerReplica]
       [La] [Lb] [hStart] [hEnd] [hStep] [seed] [options]
```

| # | Argument | Défaut | Rôle |
|---|---|---|---|
| 1 | `inputDir` | obligatoire | Répertoire contenant les fichiers d'entrée format : input_XXXXX|
| 2 | `outputDir` | obligatoire | Répertoire de sortie |
| 3 | `totalThreads` | tous les cœurs | Budget total de threads |
| 4 | `threadsPerReplica` | 4 | Threads par réplique |
| 5 | `La` | 30 | Taille du réseau le long de $\mathbf{a}$ |
| 6 | `Lb` | 30 | Taille du réseau le long de $\mathbf{b}$ |
| 7 | `hStart` | 0.0 | Champ de départ (avec `--field`) |
| 8 | `hEnd` | 1.0 | Champ d'arrêt (avec `--field`) |
| 9 | `hStep` | 0.0125 | Pas du balayage en champ (avec `--field`) |
| 10 | `seed` | 12345 | Graine maîtresse des générateurs aléatoires |

### Options booléennes

| Option | Défaut | Rôle |
|---|---|---|
| `--sqw` / `--no-sqw` | `--sqw` | Calcule $S(Q,\omega)$ à la dernière température |
| `--field` / `--no-field` | `--no-field` | Active le champ magnétique |
| `--annealing` / `--no-annealing` | `--annealing` | Exécute le recuit , `--no-annealing` exige `--ac` |
| `--adaptive-cooling` / `--no-adaptive-cooling` | géométrique 0.95 | Refroidissement adaptatif |
| `--freeze-sigma` / `--no-freeze-sigma` | figé | $\sigma$ figé en production ou adaptatif |
| `--precession-at-end` / `--no-precession-at-end` | activée | Relaxation en fin de chaque température |
| `--stagnation` / `--no-stagnation` | inactive | Détecteur de stagnation |
| `--stop-therm-on-stagnation` | inactif | Arrêt anticipé de la thermalisation |
| `--stop-prod-on-saturation` | inactif | Arrêt anticipé de la production |
| `--quiet` | inactif | Supprime la barre de progression |
| `--seed-random` | inactif | Graine tirée au hasard, affichée au démarrage |
| `--spins-per-T` | inactif | Instantané des spins à chaque température |
| `--decorrelation-log` | inactif | Journal de décorrélation par température |
| `--ac` / `--no-ac` | `--no-ac` | Susceptibilité AC |

### Options à valeur

| Option | Défaut | Rôle |
|---|---|---|
| `--overrelax N` | 2 | Sur-relaxations microcanoniques par balayage, 0 désactive |
| `--minTherm N` | 1000 | Plancher de balayages de thermalisation |
| `--maxTherm N` | 200000 | Plafond de balayages de thermalisation |
| `--targetNeff N` | 400 | $N_{eff}$ visé pendant la production |
| `--gewekeSkip F` | 0.5 | Fraction initiale écartée par le test de Geweke |
| `--initTemp T` | 10 | Température initiale du recuit |
| `--endTemp T` | 0.001 | Température finale, où le bloc $S(Q,\omega)$ est mesuré |
| `--wmax X` | historique | $w_{max}$ visé, $\Delta = \pi/X$ , 0 pour une coupe statique |
| `--window W` | 20 | Durée de la fenêtre d'intégration, $\delta\omega = 2\pi/W$ |
| `--domega X` | | Résolution visée en énergie, $W = 2\pi/X$, alternative à `--window` |
| `--nSQW N` | 40 | Nombre d'instantanés (moyenne sur une réplique) de $S(Q,\omega)$, $2 \le N $ |
| `--qz X` | 0 | Composante $z$ du vecteur de diffusion en r.l.u. de $c^*$ |
| `--field-dir x,y,z` | -1,1,0 | Direction du champ magnétique |
| `--preset P` | standard | Voir ci-dessous |
| `--memoryMultiplier K` | 3 | Seuil de mémoire thermique : thermalisation $\ge K\times\tau$ avant de déclarer deux températures décorrelées |

### Presets

| Preset | targetNeff | maxTherm | maxProd | Stagnation | Arrêts anticipés |
|---|---|---|---|---|---|
| `standard` | 400 | 200000 | 400000 | non | non |
| `screening` | 100 | 50000 | 50000 | oui | oui |
| `strict` | 1000 | 500000 | 2000000 | non | non |

Le preset strict impose de plus un $N_{eff}$ minimal de 200 par barreau et un critère de
thermalisation portant sur cinquante temps de corrélation.

Les options explicites placées après `--preset` sur la ligne de commande l'emportent.

### Susceptibilité AC

| Option | Défaut | Rôle |
|---|---|---|
| `--acReplicas N` | 8 | Répliques AC indépendantes |
| `--acH0 x` | 0.02 | Amplitude du champ oscillant |
| `--acPeriods p1,p2` | 64,128,256 | Périodes d'excitation, en balayages |
| `--acTransientPeriods N` | 5 | Périodes transitoires écartées |
| `--acMeasurePeriods N` | 30 | Périodes de mesure accumulées |
| `--acTempMax T` | 2.0 | Température de départ de l'échelle AC |
| `--acTempMin T` | 0.05 | Température d'arrêt de la grille AC |
| `--acTempStep dT` | 0.05 | Pas de la grille linéaire |
| `--acFieldDir x,y,z` | champ statique | Direction propre au champ AC |

### Codes de sortie

| Code | Signification |
|---|---|
| 0 | Succès |
| 1 | Échec d'exécution |
| 2 | Erreur d'usage |

### Exemples


Recuit complet, résolution fine du spectre :

```bash
java -cp "bin;lib/*" AppParallel inputDir outputDir 1 1 60 60 0.0 0.0 0.0 12345 --endTemp 0.01 --qz 4.67 --domega 0.05
```

Echantillonnage rapide pour 1 input :

```bash
java  -cp "bin;lib/*" AppParallel inputDir outputDir 1 1 60 60 0.0 0.0 0.0 12345 --endTemp 0.01 --qz 4.67 --preset screening --nSQW 20 --domega 0.05
```

Echantillonnage rapide pour plusieurs input :

```bash
java  -cp "bin;lib/*" AppParallel inputDir outputDir 32 1 30 30 0.0 0.0 0.0 --seed-random --endTemp 0.05 --qz 0.0 --preset screening --nSQW 20 --domega 0.05 --quiet
```


## Sorties

Chaque réplique écrit dans un sous-répertoire horodaté de `outputDir`, nommé d'après le
fichier d'entrée et le champ statique.

| Fichier | Contenu |
|---|---|
| `structure_dipolaire_h<h>.avro` | $S(Q,\omega)$ dipolaire, écrit à la dernière température du recuit |
| `structure_couleur_h<h>.avro` | Canal couleur, si terme quartique actif |
| `structure_dimere_h<h>.avro` | Canal dimère, si terme quartique actif |
| `spins_finaux_h<h>.avro` | Configuration finale de spins |
| `spins_T<T>_h<h>.avro` | Instantané de spins par température, avec `--spins-per-T` |
| `equilibration.csv` | Observables par température : énergie, aimantation, $\sigma$, $N_{eff}$, $\tau$, acceptance |
| `ac_susceptibility_<fichier>_h<h>.csv` | Résultats de la susceptibilité AC, avec `--ac` |
| `decorrelation.csv` | Journal de décorrélation par température, avec `--decorrelation-log` |


