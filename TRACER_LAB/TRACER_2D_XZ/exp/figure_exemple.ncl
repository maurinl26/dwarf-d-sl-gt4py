;####################  SCRIPT NCL POUR TRACER DES CHAMPS CALCULÉS  #####################;
;#                                                                                     #;
;# auteurs : L.Auger, F.Voitus, C.Colavolpe, B.Menetrier                               #;
;# sujet   : Permet de visualiser les résultats du dossier ${experience} via           #;
;#           la commande suivante :                                                    #;
;#                               ncl figure_hot.ncl                                    #;
;#                                                                                     #;
;#######################################################################################;

;------------  Téléchargement de données pour les graphes  ------------;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;----------------------------------------------------------------------;
begin
;------------------------  Lecture des champs  ------------------------;
NLEV  = numAsciiRow("./HOT/CRITIQUE/UBpreF/ARK/res/Z"+sprinti("%0.11i",0))
NXPT  = numAsciiCol("./HOT/CRITIQUE/UBpreF/ARK/res/Z"+sprinti("%0.11i",0))
z = asciiread("./HOT/CRITIQUE/UBpreF/ARK/res/Z"+sprinti("%0.11i",0),(/NLEV,NXPT/),"double")
x = asciiread("./HOT/CRITIQUE/UBpreF/ARK/res/X"+sprinti("%0.11i",0),(/NLEV,NXPT/),"double")
filename1 = "./HOT/CRITIQUE/UFpreF/UJ/res/T00000004000"
filename2 = "./HOT/CRITIQUE/UBpreF/ARK/res/T00000002400"
filename3 = "./HOT/CRITIQUE/UFpreF/Trap/res/T00000003750"
RES1 = asciiread(filename1,(/NLEV,NXPT/),"double")
RES2 = asciiread(filename2,(/NLEV,NXPT/),"double")
RES3 = asciiread(filename3,(/NLEV,NXPT/),"double")
plot_3  = new(3,graphic)
;----------------------------------------------------------------------;

;----------------  Valeurs des paramètres des graphes  ----------------;
map = True

; Enlève des infos inutiles
map@cnInfoLabelOn = False

map@gsnFrame = False
map@gsnDraw  = False

; Axes abscisse et ordonnée
map@cnFillMode = "AreaFill"
map@trGridType = "TriangularMesh"
map@sfXArray             = x
map@vpWidthF             = 0.6
;map@tiXAxisString       = "Distance relative (en km)"
;map@tiXAxisFontHeightF  = 0.018
map@tmXBMinorPerMajor    = 1
map@tmXBMode             = "Manual"
map@tmXBLabelFontAspectF = 1.5
map@tmXBLabelFontHeightF = .015
map@trXMinF              = -100
map@tmXBTickStartF       = -100
map@tmXBTickSpacingF     = 50
map@trXMaxF              = 200
map@tmXBTickEndF         = 200

map@sfYArray             = z
map@vpHeightF            = 0.15
;map@tiYAxisSide         = "Left"
;map@tiYAxisString       = "Altitude (en km)"
;map@tiYAxisFontHeightF  = 0.018
map@tmYLMinorPerMajor    = 1
map@tmYLMode             = "Manual"
map@tmYLLabelFontAspectF = 1.5
map@tmYLLabelFontHeightF = .015
map@trYMinF              = 0
map@tmYLTickStartF       = 1
map@tmYLTickSpacingF     = 2
map@trYMaxF              = 9.5
map@tmYLTickEndF         = 9.5

; Légende
map@lbLabelBarOn               = False
map@lbOrientation              = "vertical" 
map@cnLineLabelBackgroundColor = "transparent"
map@lbBoxLinesOn               = False
map@lbLabelAutoStride          = False
map@cnLabelBarEndStyle         = "ExcludeOuterBoxes" 
map@lbLabelFontHeightF         =.00925

; Marque les ligne de séparation sans afficher la valeur
map@cnLinesOn           = False
map@cnLineLabelPerimOn  = False
map@cnLabelMasking      = False
map@cnLevelFlags        = True
map@cnLineThicknessF    = 1
map@cnLineLabelsOn      = False
map@gsnContourNegLineDashPattern = 0
;----------------------------------------------------------------------;

;------------------------  Choix des couleurs  ------------------------;
map@cnFillOn = True

wks  = gsn_open_wks("epsi","hot")
gsn_define_colormap(wks,"temp_19lev")

map@cnLevelSelectionMode = "ExplicitLevels"
map@cnLevels     = (/ -0.003,-0.0027,-0.024,-0.0021,-0.0018,-0.0015,-0.0012,-0.0009,-0.0006,-0.0003,0.00045,0.0006,0.0009,0.0012,0.0015,0.0018,0.0021,0.0024,0.0027,0.003 /)
  map@cnFillColors = (/ 2,3,3,4,4,5,6,7,8,9,0,14,15,16,17,18,19,19,20,20,21 /)
;----------------------------------------------------------------------;

;--------------------  Paramétrisation du relief  ---------------------;
plot_3(0) = gsn_csm_contour(wks,RES1(:,:),map)
plot_3(1) = gsn_csm_contour(wks,RES2(:,:),map)
plot_3(2) = gsn_csm_contour(wks,RES3(:,:),map)

; Panel
resP = True
resP@gsnFrame = False

resP@gsnPanelLabelBar = True
  resP@lbLabelBarOn   = True
  resP@lbBoxLinesOn   = False
resP@lbOrientation  = "vertical"
  resP@pmLabelBarOrthogonalPosF = 0.025
  resP@pmLabelBarWidthF  = 0.05
  resP@pmLabelBarHeightF = 0.85

resP@cnLevelSelectionMode = "ExplicitLevels"
  resP@cnExplicitLabelBarLabelsOn = False 
  resP@cnFillOpacityF = 1.0
  resP@cnFillColors = (/ 2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,24,26,28,0,37,39,41,43,45,47,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65 /)
  resP@cnExplicitLabelBarLabelsOn = True 
  resP@lbLabelStrings             = (/"-30","-27","-24","-21","-18","-15","-12","-9","-6","-3","3","6","9","12","15","18","21","24","27","30"/)
;----------------------------------------------------------------------;

;-----------------------  Création des graphes  -----------------------;
gsn_panel(wks,plot_3,(/3,1/),resP)

frame(wks)
delete(wks)
;----------------------------------------------------------------------;

;----------------------  Visionnage des graphes  ----------------------;
system("epstopdf hot.epsi")
system("evince hot.pdf")
system("rm -f *.epsi *.pdf")
;----------------------------------------------------------------------;

end
