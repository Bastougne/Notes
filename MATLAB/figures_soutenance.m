function figures_soutenance(lots_ch5, scripts)
% FIGURES_SOUTENANCE  Les figures du manuscrit au format des diapos de soutenance.
%
%   figures_soutenance                     les cinq scripts, lots du chapitre 5 lus dans
%                                          resultats/ch5_v2/lots
%   figures_soutenance(lots_ch5)           les lots du chapitre 5 lus dans ce dossier
%   figures_soutenance(lots_ch5, scripts)  seulement ces scripts, cellule de noms ; [] pour
%                                          lots_ch5 garde le dossier par defaut
%
% Chaque script de figures du manuscrit (figures_chapitre_2, figures_chapitre_3,
% figures_campagnes, figures_chapitre_5, figures_schemas) trace aussi, quand il trouve la structure diapo
% dans son espace de travail, les versions diapo des figures que Latex/Soutenance/main.tex
% insere : memes donnees, memes calculs et meme palette, a la taille d'insertion et sans
% mise a l'echelle, donc avec un texte a diapo.fs points (7 pt pour les schemas du
% chapitre 2). Elles vont dans Latex/Soutenance/Images.
%
% LE MANUSCRIT N'EST PAS TOUCHE. Les figures du manuscrit, que chaque script trace quand
% meme, partent dans un dossier temporaire (diapo.controle), ses tableaux et
% chiffres_ch5.txt dans un autre (diapo.tab). A la fin, chaque figure de controle est
% comparee a celle d'Images/, date et identifiant du PDF mis a part : << identique >> dit
% que les scripts tracent encore le manuscrit tel quel, ce qui verifie les remaniements
% faits pour la soutenance (blocs de trace sortis en fonctions, gabarits facultatifs).
%
% LES LOTS DU CHAPITRE 5. figures_chapitre_5 retient, par etiquette, le lot le plus recent
% de son dossier ; tant que la file de calcul y ecrit, ses figures diapo s'ecarteraient de
% celles du manuscrit. Passer alors en lots_ch5 un dossier ou l'on a copie les lots que
% cite resultats/ch5_v2/tableaux/chiffres_ch5.txt.
%
% GEOMETRIE DES DIAPOS (main.tex) : beamer 4:3 a 9 pt, page de 128 x 96 mm, marges de
% 1 cm, donc \textwidth = 108 mm. diapo.textwidth la donne en points PostScript (1/72 de
% pouce, l'unite des points de MATLAB, et non le point TeX de 1/72,27 de pouce), soit
% 306,1 pt. Les figures sont dimensionnees en fractions de diapo.textwidth, celles de leur
% colonne dans main.tex : y changer une largeur veut dire la changer aussi dans le script
% qui trace la figure. MATLAB arrondit la page du PDF au point entier ; main.tex insere
% donc chaque figure a sa largeur, ce qui ne la met a l'echelle qu'a cet arrondi pres.

    here   = fileparts(mfilename('fullpath'));
    racine = fullfile(tempdir, 'figures_soutenance');
    base.img       = fullfile(fileparts(here), 'Latex', 'Soutenance', 'Images');
    base.controle  = fullfile(racine, 'manuscrit');
    base.tab       = fullfile(racine, 'tableaux');
    base.textwidth = 108 / 25.4 * 72;      % points de MATLAB, voir l'en-tete
    base.fs        = 7.5;
    if nargin >= 1 && ~isempty(lots_ch5), base.lots_ch5 = lots_ch5; end
    if nargin < 2 || isempty(scripts)
        scripts = {'figures_chapitre_2', 'figures_chapitre_3', 'figures_campagnes', ...
                   'figures_chapitre_5', 'figures_schemas'};
    end
    for d = {base.img, base.controle, base.tab}
        if ~exist(d{1}, 'dir'), mkdir(d{1}); end
    end

    for s = 1:numel(scripts)
        fprintf('\n===== %s =====\n', scripts{s});
        lancer(scripts{s}, base);
        close all;
    end
    comparer(base.controle, fullfile(fileparts(here), 'Images'));
end

function lancer(nom, diapo) %#ok<INUSD>
% Le script tourne dans l'espace de travail de cette fonction, ou il trouve diapo ; son
% clearvars -except diapo n'efface que les variables de cette fonction.
    run(nom);
end

function comparer(controle, images)
% Chaque figure de controle face a celle du manuscrit, date et identifiant mis a part.
    f = dir(fullfile(controle, '*.pdf'));
    fprintf('\n===== Figures du manuscrit retracees (%d) =====\n', numel(f));
    n_diff = 0;
    for i = 1:numel(f)
        cible = fullfile(images, f(i).name);
        if ~isfile(cible)
            etat = 'absente d''Images/';  n_diff = n_diff + 1;
        elseif strcmp(sans_horodatage(fullfile(f(i).folder, f(i).name)), sans_horodatage(cible))
            etat = 'identique';
        else
            etat = 'DIFFERENTE';  n_diff = n_diff + 1;
        end
        fprintf('  %-48s %s\n', f(i).name, etat);
    end
    fprintf('%d figure(s) a examiner.\n', n_diff);
end

function s = sans_horodatage(f)
% Repris de figures_campagnes.m.
    fid = fopen(f, 'r');
    s   = char(fread(fid, Inf, '*uint8')');
    fclose(fid);
    s = regexprep(s, '/(CreationDate|ModDate)\s*\([^)]*\)', '');
    s = regexprep(s, '/ID\s*\[[^\]]*\]', '');
end
