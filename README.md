# Conductivite_thermique

15/11/2022
Bonjour,

Le projet a été mis dans un github pour pouvoir faciliter le travail. J'ai écrit un code pour l'étude asymptotique:
https://github.com/trnhatnam/Conductivite_thermique/blob/etude_asymptotique/projet.py

J'ai affiché seulement les asymptotes à haute fréquence pour tous les demi-largeurs données pour éviter de surchager le graphique. J'ai inclus les asymptotes à basse fréquence (mais les lignes concernées sont mis en commentaire dans le code pour la raison ci-dessus).

Accès complet à la branche "étude asymptotique" : https://github.com/trnhatnam/Conductivite_thermique/tree/etude_asymptotique

-------------------------------------------------------

16/11/2022

oui pas mal un peu difficile à lire mais vous y etes, les methodes ne sont pas indispensables.

rappels:
thermal_freq = 2 * frequency  car la frequence thermqieu provient du 2nd harmonique voir calculs refs
T_depth = (np.sqrt(2*D / (2 * (math.pi) * thermal_freq)))/1e-6     # bien fonction de freq thermique

Le test sur le div/0 pourquoi pas on peut faire l`erreur en introduisant les frequences mais il faut renvoyer une liste de frequences facile à manipuler.
il faut aller une etape au dela avant de passer à la suite: stocker ces calculs dans une dataframe (pandas) puis sur votre disque en fichier .xls ou csv.
il faut faire ca maintenant car on pourra ainsi appeler les differents parametres par leur indice dans une seule dataframe et donc eviter les methodes inutiles.

elements de reflexion pour la suite:
1. il nous faut une dataframe pour stocker nos calculs en fonction de bh, ts et freq
2. on est multiniveaux et les indices seront importants car on va stocker des tableaux pour chaque bh et ts 
3. etudier la class MultiIndex dans pandas afin de faire le produit cartesien de nos entrées puis faites un reset des indices  , pourquoi ce type de produit?
rappelez vous: pour un bh pour ts je scanne en frequence , je passe à un autre bh, le meme ts je scanne en frequence et ainsi de suite jusque `a la fin de 
mes listes de parametres d`entrees.

on peut se voir en presentiel de 8.30-10.30 lundi 21 sinon il faut qu on se parle la semaine prochaine d`une maniere ou d`une autre avant mon depart.

----------------------------------------------------

28/11/2022

Bonjour, voici le code avec l'étude asymptotique : https://github.com/trnhatnam/Conductivite_thermique/blob/main/3omega_MAIN3_2022.py

Ce que j'ai modifié commence à partir de la ligne 76. Les modifications qui méritent d'être mentionnées sont :
- changement de l'étude fréquentiel (ligne 38)
- adaptation de f_u_vec en f_u_vec2 pour seulement faire l'étude asymptotique (ligne 78)
- le programme demandera les valeurs de b et ts (ligne 87 à 91) pour faciliter le travail
- ajout de la profondeur de pénétration comme seconde abscisse x (ligne 111)

Le programme va afficher l'asymptote en fonction de omega et des paramètres b et ts choisis.

Explications pour les fonctions :
- np.take(a, list_index) -> renvoie un array d'élements de a qui ont leur indice dans list_index
- ax1.twiny() -> permet de créer un second axe x qui partage le même axe y que le premier axe x
- omg2lamb(omg) -> renvoie la profondeur de pénétration en fonction de omega
- df.index[conditions] -> renvoie les indices des lignes dans df qui vérifient les conditions

Il manque le traçage des asymptotes à haute fréquence car ce que j'ai tracé est seulement l'asymptote à basse fréquence. Je vais laisser les autres prendre soin de faire cette partie si on doit la faire aussi.

TRINH Nhat-nam MAIN3

---------------------------------------------------

29/11/2022

Bonjour, oui mais pourquoi demander des entrées? il ne s`agit pas de choisir en amont mais de tracer en fonction de b et ts de telle sorte à observer les zones d`interet et faire un choix judicieux.
c est d`ailleurs l`interet d`avoir vectorisé. completer avec les autres asymptotes et tracer l`ensemble pour une valeur de ts (car toutes les ts en 2D cela encombrera le graph, on essaiera en 3D ou contour plus tard)

---------------------------

1/12/2022

voir librairie MeijerG python

    val1 = (-j*P / (4*L*k*math.pi * omega_elem)) * meijerg([[1, 3 / 2], []], [[1, 1], [0.5, 0]], j * omega_elem)  #solution approximée via fnction MeijerG on recupere reel et imaginaire
    amplitude = math.sqrt(np.real(val1) ** 2 + np.imag(val1) ** 2)
    phase = math.degrees(math.atan(np.imag(val1) / np.real(val1)))
    V3omega_asympt = 0.5 * V0 * TCR * asympt                                        # calculate thrid harmonic from DT
    V3omega = 0.5 * V0 * TCR * val1                                                 # calculate thrid harmonic from DT


- lire librairie MeijerG
- à faire : superposer la solution MeijerG sur les asymptotes --  tracer la figure
- tracer Temp vs freq elect, Temp vs freq ther, V3omega, amplitude et phase vs freq (1 axe y pour phase et 1 axe y pour amplitude) 
- faire solution integration numerique: Simpson

prevoir un Zoom le 6/12  avancer sur integration numerique

avantle 16/12 tenter de faire les courbes


-----------------------------------------

5/12/2022

test_asymptotev2 a été mis à jour (ligne 55 à fin). J'ai superposé les résultats de meijerg avec les asymptotes à basse fréquence (in phase et out phase) pour un seul b (vous pouvez changer la ligne 84 pour faire un parcours pour tous les b mais ça deviendra illisible).

J'ai retiré la profondeur de pénétration comme le b est susceptible d'évoluer, je ne vois pas comment faire pour mettre la profondeur de pénétration en second 
abscisse.

TRINH Nhat-nam MAIN3

----------------------------------------------------------------

5/12/2022
pas mal ca prend forme. tacher de comprendre lìmplementation des fonctions hypergeometriques en quoi elles sont interessantes pour nous (solution analytique de
l`integrale complexe à comparer avec Simpson integration numerique).
remarque: pourquoi dites vous que T_depth ne peut pas etre placée ? elle ne depend pas de b juste de la frequence et de la diffusivité constante intrinseque au materiau.
je propose que nous fassions le suivi de projet sur le Github vous pouvez svp transferer tous les elements dessus refs codes..etc. je mettrai mes commentaires en issue#

-------

13/12/2022
Il me semble qu'on peut seulement mettre du code sur github donc j'ai mis les deux fichiers .py
https://github.com/trnhatnam/Conductivite_thermique/tree/main
J'aurai besoin de votre nom github pour pouvoir vous ajouter dans la liste des collaborateurs.

------------

13/12/2022
voici mon login: mpbt2022
il faut que l`on fasse un point avant la fin de semaine

27/12/2022
- redaction du rapport attention sources
- partage des taches pour la redaction et la presentation :coordination importante et equilibre
- finaliser 2 courbes 2D cas extremes : une hors epaisseur substrat et une dedans
- essayer d integrer Simpsonen lieu et place de MeijerG sinon juste prendre MeijerG et Simpson et les comparer juste T vs Freq pour 2 bh et 1 ts pour observer la precision entre les 2
- pour Tfilm mince essayer de le comprendre mais on reprend ca le 3 ensemble.

---------------

04/01/2023
- La rédaction du rapport a bien avancée
- Courbes 2D étudiées en cas extrême
- Pour Simpson, comment recupérer le code de simpson_rule_MAIN3 afin de l'insérer dans f_u_vec et ainsi vectoriser le résultat à la place de MeijerG?

Rizlaine ZAROUAL MAIN3

----------------
05/01/2023
- Ajout de meijerg_only.py et simpson_only.py pour la comparaison meijerg/simpson L'objectif est de créer un csv pour meijerg et un csv pour simpson puis comparer les résultats.
