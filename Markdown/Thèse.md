## But :

- concevoir des algorithmes de navigation avec acquisition des données par magnétomètres par filtrage particulaire basés sur un formalisme SLAM

## Liens utiles :

- [[Algorithmes de Navigation par Champ Magnétique Terrestre pour des Véhicules Aérospatiaux.pdf]]
- [[Présentations]]
- documents [[à lire]]
- [[Thèses en lien]] avec mon sujet
- [[Scripts MATLAB]]

## Estimation bayésienne :

- [[Filtre optimal]]
- 3 grandes méthodes :
	- [[Filtre de Kalman]]
	- [[Filtrage particulaire]]
	- méthodes par grille (_grid-based methods_)

## Navigation inertielle :

- autonome et fiable
- mesures imparfaites (biais, facteurs d’échelle, bruits, ...)
- donc dérive croissante avec le temps
- ok pour des missions de courte durées

## Fusion de données :

- besoin de mesures externes de recalage => filtres de fusion :
	- souvent par GPS
		- complémentaires avec les centrales inertielles
		- facile à brouiller
		- pas tjs accessible
	- lidar

## Navigation magnétique :

- essor de la navigation par champ magnétique en robotique terrestre (indoor, piéton)
- de + en + envisagée pour l'aéroporté car :
	- amélioration des magnétomètres
	- amélioration des calibrations de magnétomètres en milieu perturbé magnétiquement
	- indépendant de la météo
	- possible au dessus de la terre et de la mer
	- capteurs compacts et embarquables

## Exploitation du champ :

- 3 manières d'exploiter le champ magnétique terrestre :
	- comparaison d'une mesure locale avec une carte spatiale préexistante précise
		- déduction des paramètres cinématiques du véhicule : position, vitesse, altitude, ...
		- le plus connu dans la littérature
	- uniquement avec des mesures locales :
		- obtention de certains paramètres cinématiques : vitesses rectiligne et angulaire
	- construction d'une carte au voisinage de la trajectoire en cours de mission à partir d'un modèle peu ou mal connu :
		- assimilé aux SLAM (initialement par vision artificielle)

## Fonction d'observation :
- idéalement, on embarque une carte qui modélise $h_k$
- en général, on ne la connait pas, ou pas assez précisément (surtout en navigation magnétique)
- on ne peut pas calculer les innovations successives des filtres
- on ne peut donc pas calculer le gain de Kalman (filtres de Kalman), ni mettre à jour les poids des particules (filtres particulaires)
- on va donc interpoler des données connues d'avance (potentiellement bruitées) par [[Krigeage ordinaire]], et utiliser cette interpolation comme fonction d'observation

## SLAM magnétique :

- déjà présents dans la littérature
- mais limités à des variabilités importantes du champ magnétique
- navigation aéroportée 100% magnétique très peu présente
- [[SLAM]]
