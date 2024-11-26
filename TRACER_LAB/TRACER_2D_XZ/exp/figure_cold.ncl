;####################  SCRIPT NCL POUR TRACER DES CHAMPS CALCULÉS  #####################;
;#                                                                                     #;
;# auteurs : L.Auger, F.Voitus, C.Colavolpe, B.Menetrier                               #;
;# sujet   : Permet de visualiser les résultats du dossier ${experience} via           #;
;#           la commande suivante :                                                    #;
;#                            ncl 'var="X"' 'step=n' figure.ncl                        #;
;#           où ${X} une des variables {U,W,T}                                         #;
;#                                                                                     #;
;#######################################################################################;

;---------------------  Definition des varialbes  ---------------------;
;experience = "COLD_QE_50"
;experience = "COLD_EE_50"
;experience = "COLD"
experience = "TEST"
dossier = "./"+experience+"/res/"
;dossier = "../res/"
;----------------------------------------------------------------------;

;------------  Téléchargement de données pour les graphes  ------------;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;----------------------------------------------------------------------;
begin
;------------------------  Lecture des champs  ------------------------;
NLEV  = numAsciiRow(dossier+"Z"+sprinti("%0.11i",0))
NXPT  = numAsciiCol(dossier+"Z"+sprinti("%0.11i",0))
z = asciiread(dossier+"Z"+sprinti("%0.11i",0),(/NLEV,NXPT/),"double")
x = asciiread(dossier+"X"+sprinti("%0.11i",0),(/NLEV,NXPT/),"double")
filename = dossier+var+sprinti("%0.11i",step)
RES = asciiread(filename,(/NLEV,NXPT/),"double")
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
map@tiXAxisString      = "x(km)"
map@tiXAxisFontHeightF = 0.01
map@vpWidthF           = 0.5
map@tmXBMinorPerMajor = 1
map@tmXBMode          = "Manual"
map@tmXBTickEndF      = 20
map@trXMaxF           = 19
map@tmXBTickSpacingF  = 2.5
;map@tmXBTickStartF    = 10000
map@trXMinF           = 0
map@tmXBLabelFontAspectF = 1.5
map@tmXBLabelFontHeightF = .0065
map@tmXBMajorOutwardLengthF = 0.0
map@tmXBMinorOutwardLengthF = 0
map@tmXBMajorLengthF = 0.0025
map@tmXBMajorThicknessF = 0.5
map@tmXBMinorLengthF = 0

map@sfYArray = z
;map@tiYAxisSide        = "Left"
map@tiYAxisString      = "z(km)"
map@tiYAxisFontHeightF = 0.01
map@vpHeightF          = 0.15
map@tmYLMinorPerMajor = 1
map@tmYLMode          = "Manual"
;map@tmYLTickEndF      = 20
map@trYMaxF           = 5.
map@tmYLTickEndF      = 5.
map@tmYLTickSpacingF  = 1
map@tmYLTickStartF    = 0
map@tmYLLabelFontAspectF = 1.5
map@tmYLLabelFontHeightF = .0065
map@tmYLMajorOutwardLengthF = 0.0
map@tmYLMinorOutwardLengthF = 0.0
map@tmYLMajorLengthF = 0.0025
map@tmYLMajorThicknessF = 0.5
map@tmYLMinorLengthF = 0


map@tmBorderThicknessF = 0.8

; Légende
map@lbLabelBarOn               = True
map@lbOrientation              = "vertical" 
map@cnLineLabelBackgroundColor = "transparent"
map@lbBoxLinesOn               = False
map@lbLabelAutoStride          = False
map@cnLabelBarEndStyle         = "ExcludeOuterBoxes"
map@lbLabelFontHeightF         =.00925
;map@lbTitleOn                  = True
;map@lbTitleString              = "x10~S~-3"
;map@lbTitleFontHeightF         = 0.025
;map@pmLabelBarOrthogonalPosF   = .001 
;map@pmLabelBarParallelPosF     = 0.51
;map@pmLabelBarWidthF   = .04

map@pmLabelBarWidthF           = .025
map@pmLabelBarHeightF          = .165
map@lbLabelFontThicknessF      = 0.8
map@lbLabelFontHeightF         =.0055

; Marque les ligne de séparation sans afficher la valeur
map@cnLinesOn           = True
map@cnLineLabelPerimOn  = False
map@cnLabelMasking      = False
map@cnLevelFlags        = False
map@cnLineThicknessF    = 0.5
map@cnLineLabelsOn      = False
map@gsnContourNegLineDashPattern = 0
;----------------------------------------------------------------------;

;------------------------  Choix des couleurs  ------------------------;
map@cnFillOn = True
wks  = gsn_open_wks("epsi",var+sprinti("%0.11i",step))

; Général de [-max ; +max]
;gsn_define_colormap(wks,"temp_19lev")
;map@cnLevelSelectionMode = "ManualLevels"
;  map@cnMaxLevelValF  = doubletointeger(10*max((/abs(max(RES(:,:))),abs(min(RES(:,:)))/)))/10.
;  map@cnMinLevelValF  = -map@cnMaxLevelValF
;  map@cnLevelSpacingF = (map@cnMaxLevelValF-map@cnMinLevelValF)/15
;  map@cnFillColors = (/ 2,3,4,5,6,7,8,9,0,14,15,16,17,18,19,20,21 /)


; COLD BUBLLE
;gsn_define_colormap(wks,"MPL_ocean")
;map@cnLevelSelectionMode = "ExplicitLevels"
;  map@cnLevels     =  (/-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1/)
;  map@cnFillColors = (/ 40,46,52,58,64,70,76,82,88,94,100,110,0 /)

gsn_define_colormap(wks,"MPL_Greys")
map@cnLevelSelectionMode = "ExplicitLevels"
  map@cnLevels     =  (/-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0.085/)
  map@cnFillColors = (/ 128,109,100,91,82,73,64,55,46,37,28,19,10,0 /)
  map@lbLabelStrings = (/ "-13","-12","-11","-10","-9","-8","-7","-6","-5","-4","-3","-2","-1","0" /)
;----------------------------------------------------------------------;

;--------------------  Paramétrisation du relief  ---------------------;
plot = gsn_csm_contour(wks,RES(:,:),map)


;----------------------------------------------------------------------;

;----------------------------------------------------------------------;

;-----------------------  Création des graphes  -----------------------;
draw(plot)
frame(wks)
delete(wks)
;----------------------------------------------------------------------;

;----------------------  Visionnage des graphes  ----------------------;
system("epstopdf "+var+sprinti("%0.11i",step)+".epsi")
system("evince "+var+sprinti("%0.11i",step)+".pdf")
system("rm -f *.epsi *.pdf")
system("ls -t "+dossier+var+"*")
;----------------------------------------------------------------------;

end
