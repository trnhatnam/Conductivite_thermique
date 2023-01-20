# Projet multidisciplinaire Polytech 3 - MAIN (2022-2023) : Mesure de la conductivité thermique de matériaux submicroniques avec la méthode des 3-omega

Libraries utilisées : numpy, matplotlib, pandas, mpmath (normalement ces libraries sont incluses dans anaconda)

Les codes importants sont :
- 3omega_plots_MAIN.py (dans le dossier Code) : vous pouvez changer les paramètres du substrat et de l'installation en changeant les variables au début du programme. Il suffit juste de lancer le code. Des csv vont être crée ce qui est normal.

- thinfilmsSiN100nm.py et thinfilmsSiN50nm.py : les paramètres ne doivent pas être changés car le programme extrait les données d'une expérience pratique dans un fichier tableau. Il suffit que les fichiers SiN100nm.xlsx et SiN50nm.xlsx soient dans le même dossier que les codes pour que les codes s'éxecutent correctement.

- comparaison.py, meijerg_only.py, simpson_only.py (dans le dossier Comparaison_meijerg_simpson) :
1/ lancer meijerg_only.py et simpson_only.py avec les mêmes paramètres (cf variables) pour créer les csv (à noter que vous pouvez changer les paramètres du modèle mais pour la comparaison il faut mettre les mêmes paramètres);
2/ une fois les csv crées vous pouvez lancer comparaison.py
