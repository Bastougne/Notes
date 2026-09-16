# Localité, multimodalité et ce que le CA-OK gagne vraiment

Note d'analyse, écrite le 2026-09-15. Elle part d'une question de Bastien — *la réunion marche aussi bien que le clustering ; si le champ était vraiment non stationnaire, le CA-OK, qui ajuste ses hyperparamètres par mode, n'aurait-il pas dû gagner ?* — et y répond avec les deux lots à 1000 essais qui portent les trois méthodes. Elle sert trois choses : garder les chiffres, **justifier plusieurs choix du chapitre**, et poser les expériences qui rendraient multimodalité et instationnarité visibles.

Les chiffres sont mesurés et citables. Les **hypothèses** sont signalées comme telles : aucune n'est testée.

## Les lots et les métriques

- `resultats/20260912_030842_reference_libre_1000runs.mat` — les quatre modèles contre la maille, scénario `reference`, $\sigma_0 = 8$ km, six mailles de 2,5 à 7,5 km.
- `resultats/20260914_235935_reference_libre_1000runs.mat` — le régime multimodal, $\sigma_0 = 12$ km, cinq mailles de 4,5 à 7,5 km, avec le statique et la carte connue.

RPF, 5000 particules, 150 pas, $\alpha_{\mathrm{cov}} = 2$, $\tilde{N}_{\min} = 25$, $\tilde{N}_{\max} = 2000$, bande passante du mean-shift 1000 m. **acq** : erreur médiane au pas 30. **fin** : médiane sur les essais de la médiane des dix derniers pas. **Mahal.** : taux de convergence au critère de Mahalanobis (99 %). **multi** : part des couples (pas, essai) où le mean-shift du CA-OK rend plus d'un mode. **pts** : nombre moyen de points retenus. Script de dépouillement : `scratchpad/analyse_cak_reunion.m` de la session du 2026-09-15 (temporaire ; il tient en quarante lignes et se réécrit).

## Les chiffres

### Scénario de référence, $\sigma_0 = 8$ km

| maille | AK acq / fin / Mahal. | CA-OK acq / fin / Mahal. (multi) | réunion acq / fin / Mahal. (pts) |
|---|---|---|---|
| 2,5 km | 364 / 239 / 100 % | 365 / 242 / 100 % (7 %) | 386 / 307 / 100 % (22) |
| 3,5 km | 473 / 381 / 100 % | 475 / 379 / 100 % (10 %) | 590 / 488 / 99 % (18) |
| 4,5 km | 598 / 663 / 99 % | 608 / 663 / 99 % (16 %) | 748 / 757 / 97 % (15) |
| 5,5 km | 772 / 1035 / 90 % | 812 / 1052 / 91 % (26 %) | 828 / 1086 / 86 % (14) |
| 6,5 km | 861 / 1357 / 87 % | 942 / 1388 / 87 % (39 %) | 926 / 1330 / 81 % (14) |
| 7,5 km | 1459 / 1684 / 87 % | 1699 / 1720 / 87 % (54 %) | 1450 / 1813 / 74 % (13) |

Au critère de seuil (2000 m), les trois restent à trois points les uns des autres à toutes les mailles.

### Régime multimodal, $\sigma_0 = 12$ km

| maille | statique acq / fin / Mahal. | AK | CA-OK (multi) | réunion (pts) |
|---|---|---|---|---|
| 4,5 km | 682 / 773 / 99 % | 627 / 667 / 98 % | 640 / 664 / 98 % (17 %) | 828 / 746 / 97 % (22) |
| 5,0 km | 692 / 910 / 95 % | 643 / 802 / 94 % | 665 / 804 / 94 % (23 %) | 784 / 803 / 92 % (21) |
| 5,5 km | 816 / 1142 / 89 % | 881 / 1042 / 90 % | 906 / 1050 / 90 % (28 %) | 944 / 1087 / 86 % (20) |
| 6,5 km | 851 / 1474 / 83 % | 1003 / 1375 / 86 % | 1162 / 1397 / 86 % (42 %) | 1055 / 1416 / 79 % (19) |
| 7,5 km | 1297 / 1704 / 78 % | 1673 / 1726 / 86 % | 2129 / 1790 / 84 % (57 %) | 1679 / 1808 / 75 % (18) |

## Ce que disent les chiffres

1. **Le CA-OK n'est meilleur que l'AK nulle part.** Aux mailles lâches il est **pire à l'acquisition** — +9 et +16 % à 6,5 et 7,5 km en référence, +16 et +27 % en multimodal — et rejoint l'AK en fin de vol. C'est précisément là qu'il est multimodal le plus souvent (39 à 57 % des pas).
2. **La réunion n'est pas aussi bonne.** Aux mailles fines, 14 à 28 % d'erreur en plus (307 contre 239 m à 2,5 km, 488 contre 381 à 3,5 km). Aux mailles lâches, elle retrouve la précision de l'AK mais perd jusqu'à **13 points de Mahalanobis** (74 contre 87 % à 7,5 km).
3. **L'AK fait déjà ce qu'on reproche à la réunion, sans en souffrir.** Pendant l'acquisition sa fenêtre unique, de rayon $6\sqrt{\lambda_{\max}}$ — 48 km au départ avec $\sigma_0 = 8$ km —, couvre tous les modes avec **un seul ajustement et une seule moyenne**. Il ne perd rien contre le CA-OK : les hyperparamètres varient donc peu à l'échelle qui sépare les modes, de l'ordre de $\sigma_0$.
4. **Le champ est non stationnaire, mais à grande échelle.** En multimodal, l'AK bat le statique de **7 à 14 % en erreur finale** de 4,5 à 6,5 km, sa NEES finale tombe de 0,67 à 0,42 à 4,5 km, et il gagne 8 points de Mahalanobis à 7,5 km. Les hyperparamètres locaux servent ; leur localité *par mode* n'ajoute rien.
5. **Mais le statique acquiert mieux que l'AK dès 5,5 km en multimodal** (816 contre 881 m, 851 contre 1003, 1297 contre 1673), avant de finir moins bien. Pendant l'acquisition, le fenêtrage coûte quelque chose que la fin de vol rembourse.

## Pourquoi, d'après le code

- **Sur un pas unimodal, le CA-OK est l'AK** : `query_cak` appelle `query_ak` quand le mean-shift rend un mode. C'est 90 % des pas aux mailles fines ; le CA-OK ne peut s'écarter de l'AK que sur les pas multimodaux.
- **La réunion** (`query_boules`) prend les quatre plus proches échantillons de chaque particule, sans plancher, et fait **un seul** ajustement et **une seule** résolution sur leur union — 13 à 22 points en moyenne, sous $\tilde{N}_{\min}$.
- **Là où le CA-OK est multimodal, il ne peut pas être local.** Une fenêtre au plancher de 25 points couvre un disque de rayon $\approx 2{,}8\,\Delta$ : 7 km à 2,5 km de maille, **21 km à 7,5 km**, soit quatre longueurs de corrélation ($\ell \approx 5$ km). La séparation $\ell \ll r \ll L_\theta$ de la note de stationnarité est intenable aux mailles lâches — qui sont justement celles où le CA-OK se partitionne.

## Hypothèses, non testées

- **H1 — la réunion perd la cohérence par son nombre de points** (critère 2), et non parce qu'elle enjambe les modes (critère 3). Point 3 à l'appui : l'AK enjambe les modes lui aussi, mais avec au moins 25 points, et ne perd rien.
- **H2 — le CA-OK perd à l'acquisition parce que le mean-shift fragmente des nuages larges mais unimodaux.** 54 % de pas « multimodaux » à 7,5 km dès $\sigma_0 = 8$ km, c'est beaucoup pour un vrai multimodal. Des morceaux d'un même mode sont alors krigés sous des ajustements différents, sur des fenêtres plus petites.
- **H3 — réajuster les hyperparamètres à chaque pas injecte un bruit d'estimation qui coûte pendant l'acquisition**, et la localité ne paie qu'une fois le nuage resserré. Point 5 à l'appui, et l'ordre statique < AK < CA-OK des erreurs d'acquisition aux mailles lâches suit l'ordre des effectifs d'ajustement.
- **H4 — les hyperparamètres varient à une échelle $L_\theta$ bien plus grande que l'écart entre modes**, mais plus petite que la carte.

## Ce que ça justifie dans le chapitre

- **$\tilde{N}_{\min}$** : si H1 tient, la réunion est la démonstration empirique du critère 2 — elle perd jusqu'à 28 % d'erreur et 13 points de cohérence pour être passée sous le plancher. C'était déjà son rôle déclaré dans `notes_campagne_E.md` (« le second échec motive $n_{\min}$ »).
- **Le réajustement local** (critère 3) : l'AK bat le statique en fin de vol en multimodal. La localité se justifie à l'échelle de dizaines de kilomètres, pas à celle d'un mode.
- **Le CA-OK se justifie par le coût et par la taille de fenêtre pendant l'acquisition** (critère 4), pas par la précision ni par la localité des hyperparamètres — cohérent avec le choix de Bastien d'introduire le CAK par le plafond de fenêtre. **Présenter la réunion comme la preuve qu'il faut clusteriser serait fragile.**
- **La double inégalité $\ell \ll r \ll L_\theta$** reçoit une conséquence mesurable : aux mailles lâches, la fenêtre minimale est déjà trop grande pour être locale.
- **Le critère de Mahalanobis** est à nouveau celui qui départage les méthodes ; le seuil ne voit presque rien.

## Rendre la multimodalité et l'instationnarité visibles

Par ordre de coût, du moins cher au plus démonstratif. Les commandes suivent l'idiome de `notes_campagne_E.md` et **n'ont pas été lancées** ; elles passent après les lots prioritaires de la file.

### A. Mesurer d'abord, sans filtre

1. **Carte des hyperparamètres locaux.** Ajuster $\ell$ et $\sigma_f$ par maximum de vraisemblance sur des fenêtres de 25, 50 et 250 points centrées tous les quelques kilomètres le long du corridor, puis sur toute la carte. Répéter sur plusieurs tirages de relevé pour séparer le bruit d'estimation de la variation réelle, et tracer le variogramme de $\log\ell$ : sa portée **est** $L_\theta$. Tranche H4, et dit si le corridor du scénario est simplement trop homogène. Script à écrire, une cinquantaine de lignes, quelques minutes de calcul.
2. **Distances entre modes pendant l'acquisition.** Enregistrer les centres du mean-shift sur quelques essais pour comparer l'écart entre modes à $L_\theta$ et au rayon des fenêtres.

### B. Isoler les mécanismes sur le scénario existant

1. **Hyperparamètres globaux contre locaux, AK et CA-OK** (tranche H3 et la part de la localité) — les résultats `fenetre` existent déjà dans le lot multimodal, sur les mêmes graines :
   ```bash
   matlab -batch "cd('D:/bhubert/Desktop/Notes/MATLAB'); over.campagne='libre'; over.compare.models={'ak','cak'}; over.compare.filters={'RPF'}; over.sigma_0=[12000 12000 2 2]; over.krig.refit='global'; over.balayage=struct('champ','krig.pas','valeurs',[4500 5500 7500]); over.mc.n_runs=1000; over.krig.n_surveys=100; campagne_chapitre_4"
   ```
2. **La réunion au-dessus du plancher** (tranche H1) : si le Mahalanobis remonte vers 87 % quand l'union dépasse 25 points, c'était le nombre.
   ```bash
   matlab -batch "cd('D:/bhubert/Desktop/Notes/MATLAB'); over.campagne='libre'; over.compare.models={'ak_boules'}; over.compare.filters={'RPF'}; over.krig.pas=7500; over.balayage=struct('champ','krig.k_voisins','valeurs',[4 8 16 32]); over.mc.n_runs=1000; over.krig.n_surveys=100; campagne_chapitre_4"
   ```
3. **La fragmentation du mean-shift** (tranche H2) : si l'erreur d'acquisition du CA-OK rejoint l'AK quand la bande passante grandit, c'était la fragmentation.
   ```bash
   matlab -batch "cd('D:/bhubert/Desktop/Notes/MATLAB'); over.campagne='libre'; over.compare.models={'cak'}; over.compare.filters={'RPF'}; over.krig.pas=7500; over.balayage=struct('champ','krig.bandwidth','valeurs',[1000 2000 4000 8000]); over.mc.n_runs=1000; over.krig.n_surveys=100; campagne_chapitre_4"
   ```

### C. Faire durer la multimodalité

Le premier virage lève l'ambiguïté dès le pas 33 : le multimodal ne dure pas assez pour que le CA-OK puisse s'exprimer. Deux leviers, qui existent déjà :

1. **Le scénario `acquisition`** — vol en ligne, $\sigma_0 = 15$ km, 60 pas — n'a pas de virage pour lever l'ambiguïté. À vérifier qu'il tourne encore avec le code actuel.
   ```bash
   matlab -batch "cd('D:/bhubert/Desktop/Notes/MATLAB'); over.scenario='acquisition'; over.campagne='libre'; over.compare.models={'ok','ak','cak','ak_boules'}; over.compare.filters={'RPF'}; over.mc.n_runs=1000; over.krig.n_surveys=100; campagne_chapitre_4"
   ```
2. **Un $\sigma_0$ plus grand**, 20 km, sur la tondeuse — attention à la marge au bord, commentée dans le scénario `reference`.

Et **mesurer sur la phase multimodale seule** : erreur pendant les pas multimodaux, pas de levée d'ambiguïté, part des essais qui se verrouillent sur un mauvais mode. Les métriques de fin de vol la diluent.

### D. L'expérience démonstrative : une instationnarité contrôlée

La seule qui rendrait la localité **visible à coup sûr** : un champ synthétique dont on connaît les hyperparamètres en tout point. Par exemple deux tirages Matérn 3/2 stationnaires très différents — $\ell$ de 2 et 10 km, $\sigma_f$ de 50 et 250 nT — fondus le long d'une frontière par un masque lisse, ou un noyau non stationnaire de Paciorek–Schervish avec $\ell(x)$ variant sur une échelle $L_\theta$ réglable. Placer l'incertitude initiale à cheval sur la frontière, pour que les modes tombent dans des régimes différents, puis balayer $L_\theta$ contre l'écart entre modes et le rayon des fenêtres.

Prédiction : statique < AK < CA-OK quand $L_\theta$ est de l'ordre de l'écart entre modes ; les trois confondus quand $L_\theta$ le dépasse de beaucoup. C'est la validation directe de $\ell \ll r \ll L_\theta$. Il faut une carte au format de `carte_magnetometrie_anomalie_mexique.mat` (champs `mesures` et `pas`) et pointer `par.map.file` dessus. Variante sans synthèse : choisir un corridor réel qui traverse des régions contrastées, repérées grâce à la carte de A1.

### E. Figures pour le manuscrit

- La carte de $\ell(x)$ et $\sigma_f(x)$ avec la trajectoire par-dessus : l'instationnarité d'un coup d'œil.
- Un instantané du nuage à un pas multimodal, avec le disque unique de l'AK et les disques par mode du CA-OK sur la carte : pourquoi la fenêtre de l'AK enfle, ce que fait le CA-OK. Demande d'enregistrer les particules d'un essai.
- Les performances contre $L_\theta$, si D tourne.

## Points de vigilance

- **Appariement** : les comparaisons de B supposent les mêmes graines et les mêmes tirages de relevé que les lots existants — ne toucher ni à `run.seed` ni à `n_surveys`.
- **« multi » n'est pas « multimodal »** : c'est ce que rend le mean-shift à 1000 m de bande passante. Tant que H2 n'est pas tranchée, ne pas l'appeler part de pas multimodaux dans le manuscrit.
- **Noms des modèles** : `OK adaptatif`, `CA-OK`, `OK par réunion`, `OK éclairci`, `OK carte totale`, `carte connue`. Toujours apparier sur des motifs exclusifs.
- **Le temps de calcul n'est pas une contrainte** : 1000 essais d'emblée.
