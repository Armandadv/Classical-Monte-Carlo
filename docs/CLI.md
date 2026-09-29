# Référence CLI — `AppParallel`

> Page de la documentation — [retour au README](../README.md)

Monte Carlo classique de spins (Kitaev-Heisenberg anisotrope + Zeeman + terme quartique
optionnel), parallélisé en damier, dynamique Landau-Lifshitz pour S(Q,ω), sorties Avro.

```
java -cp "bin;lib/*" AppParallel <inputDir> <outputDir> [totalThreads] [threadsPerReplica]
       [La] [Lb] [hStart] [hEnd] [hStep] [seed] [options]

java -cp "bin;lib/*" AppParallel --bench <La> <Lb> [J1 J2 J3]      # benchmark synthétique
java -cp "bin;lib/*" AppParallel --bench <fichierInput> <La> <Lb>  # benchmark modèle réel
java -cp "bin;lib/*" AppParallel -h | --help
```

> Windows : séparateur de classpath `;`. Linux/macOS : `:`. Sur machine partagée, épingler
> la JVM : `-XX:ActiveProcessorCount=N`.

## Arguments positionnels

| # | Argument | Défaut | Rôle |
|---|---|---|---|
| 1 | `inputDir` | — (obligatoire) | Répertoire contenant les fichiers d'entrée `bcaoExplor_*` |
| 2 | `outputDir` | — (obligatoire) | Répertoire de sortie (créé s'il faut) |
| 3 | `totalThreads` | tous les cœurs | Budget total de threads (`concurrentReplicas = totalThreads / threadsPerReplica`) |
| 4 | `threadsPerReplica` | `4` | Threads par réplique (choisir avec `--bench`) |
| 5 | `La` | `30` | Taille du réseau le long de **a** |
| 6 | `Lb` | `30` | Taille du réseau le long de **b** |
| 7 | `hStart` | `0.0` | Champ de départ (uniquement avec `--field`) |
| 8 | `hEnd` | `1.0` | Champ d'arrêt (uniquement avec `--field`) |
| 9 | `hStep` | `0.0125` | Pas du balayage en champ (uniquement avec `--field`) |
| 10 | `seed` | `12345` | Graine maîtresse (rejoue exactement le même run) |

Sans `--field`, les trois arguments de champ ne sont ni lus ni validés : une seule réplique à h = 0.

## Options booléennes (paires `--x` / `--no-x`)

| Option | Défaut | Rôle |
|---|---|---|
| `--sqw` / `--no-sqw` | `--sqw` | Calcule S(Q,ω) à la dernière température du recuit |
| `--field` / `--no-field` | `--no-field` | Active la boucle en champ `hStart → hEnd` par pas `hStep` |
| `--annealing` / `--no-annealing` | `--annealing` | Exécute le recuit ; `--no-annealing` exige `--ac` (`EXIT_USAGE` sinon) |
| `--adaptive-cooling` / `--no-adaptive-cooling` (alias : `--cooling`) | géométrique ×0,95 | Refroidissement adaptatif de Huang (pas ∝ σ_E, resserré près de la transition) |
| `--freeze-sigma` / `--no-freeze-sigma` (alias du `--no-` : `--legacy-sigma`) | `--freeze-sigma` | σ figé pendant la production (chaîne markovienne) ou adaptatif (historique) |
| `--precession-at-end` / `--no-precession-at-end` | `--precession-at-end` | Relaxation Landau-Lifshitz (`Precession1D`) en fin de chaque température |
| `--stagnation` / `--no-stagnation` | `--no-stagnation` | Détecteur de stagnation (diagnostic ; arêts anticipés à activer séparément) |
| `--stop-therm-on-stagnation` / `--no-…` | désactivé | Coupe la thermalisation quand le détecteur conclut à la stagnation (exige `--stagnation` ou `--preset screening`) |
| `--stop-prod-on-saturation` / `--no-…` | désactivé | Coupe l'extension de production quand le gain de N_eff sature (même exigence) |
| `--quiet` | désactivé | Supprime barre de progression et résumé par température |
| `--seed-random` | désactivé | Graine tirée au hasard (affichée au démarrage ; la redonner en positionnel 10 reproduit le run) |
| `--spins-per-T` | désactivé | Instantané des spins (Avro) après la production de **chaque** température |
| `--decorrelation-log` | désactivé | Journal de décorrélation par température (τ, g, N_eff) |
| `--ac` / `--no-ac` | `--no-ac` | Susceptibilité AC χ′/χ″ en plus du recuit (ou à sa place avec `--no-annealing`) |

## Options à valeur

| Option | Défaut | Rôle |
|---|---|---|
| `--overrelax N` | `2` | N sur-relaxations microcanoniques par balayage Metropolis (`0` = désactivé ; interdit avec une ligne `QUARTIC` dans l'input) |
| `--minTherm N` | contrôleur (`1000`) | Plancher de sweeps de thermalisation (`minSweeps`) |
| `--maxTherm N` | contrôleur (`200000`) | Plafond de sweeps de thermalisation (`maxSweeps`) |
| `--targetNeff N` | `400` | N_eff visé pendant la production (étend les sweeps jusqu'à l'atteindre, dans la limite du plafond) |
| `--gewekeSkip F` | `0.5` | Fraction initiale écartée par le test de Geweke |
| `--memoryMultiplier X` | — | Multiplicateur de mémoire des tampons d'équilibration |
| `--initTemp T` | `10` | Température initiale du recuit |
| `--endTemp T` | `0.001` | Température finale ; **le bloc S(Q,ω) est mesuré à cette température**. Exige `0 < endTemp < initTemp` |
| `--wmax X` | historique | w_max visé (unités de Spin) → Δ = π/X (Nyquist) ; `X = 0` : mode statique (coupe unique à ω = 0) |
| `--window W` | `20` | Durée de la fenêtre d'intégration par instantané S(Q,ω) ; fixe δω = 2π/W |
| `--domega X` | — | Résolution visée en énergie (ħ = 1) : **écriture alternative** de `--window` (W = 2π/X). Incompatible avec `--window` (`EXIT_USAGE` si fournis ensemble) ; avec `--wmax > 0`, exige `X < wmax` |
| `--nSQW N` | `40` (historique) | Nombre d'instantanés du bloc S(Q,ω) (N+1 mesures spectrales moyennées, Avro écrit à la (N+1)-ième). Exige `2 ≤ N ≤ 1 000 000` ; coût du bloc proportionnel à N |
| `--qz X` | `0` | Composante z du vecteur de diffusion (r.l.u. de c*) ; `0` = réseau 2D strict (pairing miroir actif) ; valeur BaCoAsO historique : `4.67` (désactive le pairing) ; tracée dans la propriété `Qzz` des Avro de structure |
| `--field-dir x,y,z` | `-1,1,0` | Direction du champ dans le tenseur g (vecteur non nul) |
| `--preset P` | `standard` | Voir « Presets » ci-dessous |

### Notes S(Q,ω)

- La largeur de raie élastique est **automatique** : FWHM = 1 bin quelle que soit la
  résolution (fenêtre gaussienne σ = 0.375·ndt — plus rien à régler).
- La largeur visible après convolution instrumentale (Lorentzienne γ, p.ex. CAMEA) suit le
  profil de Voigt : `F ≈ 0.5346·2γ + √(0.2166·(2γ)² + δω²)` — plancher instrumental 2γ.
- Les canaux `structure_couleur` / `structure_dimere` ne sont calculés et écrits **que si
  le terme quartique est actif** (ligne `QUARTIC` dans l'input, b ≠ 0).

## Susceptibilité AC (`--ac`)

| Option | Défaut | Rôle |
|---|---|---|
| `--acReplicas N` | `8` | Répliques AC indépendantes par (fichier, champ statique) |
| `--acH0 x` | `0.02` | Amplitude du champ oscillant (x > 0) |
| `--acPeriods p1,p2,…` | `64,128,256` | Périodes d'excitation, en sweeps (entiers ≥ 2, croissantes) |
| `--acTransientPeriods N` | `5` | Périodes de régime transitoire écartées avant la démodulation |
| `--acMeasurePeriods N` | `30` | Périodes de mesure accumulées dans la démodulation |
| `--acTempMax T` | `2.0` | Température de départ (la plus haute) de chaque échelle de refroidissement AC |
| `--acTempMin T` | `0.05` | Température d'arrêt (la plus basse) de la grille AC |
| `--acTempStep dT` | `0.05` | Pas de la grille linéaire de température AC |
| `--acFieldDir x,y,z` | `--field-dir` | Direction propre au champ AC (défaut : celle du champ statique) |

Chaque période de `--acPeriods` rejoue une échelle de refroidissement **complète** depuis
`--acTempMax` — le coût total est proportionnel à `Σᵢ pᵢ × (échelle)`. Résultats dans
`equilibration.csv` et fichiers AC dédiés.

## Presets

| Preset | targetNeff | maxTherm | maxProd | Stagnation | Arrêts anticipés |
|---|---|---|---|---|---|
| `standard` (défaut) | 400 | 200 000 | 400 000 | non | non |
| `screening` | 100 | 50 000 | 50 000 | **oui** | **oui** (therm + prod) |
| `strict` | 1000 | 500 000 | 2 000 000 | non | non (minNeff 200, τ×50) |

Les options explicites placées **après** `--preset` sur la ligne de commande l'emportent
(parseur linéaire, « dernier lu gagne »).

## Codes de sortie

| Code | Constante | Signification |
|---|---|---|
| `0` | `EXIT_OK` | Succès (aussi `-h`/`--help`) |
| `1` | `EXIT_FAILURE` | Échec d'exécution (I/O, run interrompu…) |
| `2` | `EXIT_USAGE` | Erreur d'usage : argument invalide, conflit d'options, chemins absents |

## Exemples

```bash
# Recuit complet, modèle réel, 60×60, champ nul, résolution fine des Bragg
java -cp "bin;lib/*" AppParallel benchmark/modele_1/input sortie 6 6 60 60 0.0 0.0 0.0 12345 \
    --endTemp 0.01 --qz 4.67 --domega 0.05

# Criblage rapide (1 cœur) : N_eff 100, arrêts anticipés, bloc spectral à 20 instantanés
java -XX:ActiveProcessorCount=1 -cp "bin;lib/*" AppParallel in out 1 1 60 60 0.0 0.0 0.0 12345 \
    --endTemp 0.01 --qz 4.67 --preset screening --nSQW 20

# Susceptibilité AC seule (sans recuit)
java -cp "bin;lib/*" AppParallel in out 6 6 12 12 0.0 0.0 0.0 7 --no-annealing --ac \
    --acH0 0.05 --acPeriods 128,512,2048,8192 --acMeasurePeriods 480

# Choisir threadsPerReplica sur un réseau donné
java -cp "bin;lib/*" AppParallel --bench benchmark/modele_1/input 60 60
```
