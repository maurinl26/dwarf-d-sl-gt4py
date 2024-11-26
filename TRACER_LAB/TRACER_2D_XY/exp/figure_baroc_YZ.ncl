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
;experience = "EE_baroc"
experience = "TEST"
dossier = "./"+experience+"/res/X/"
;dossier = "../res/"
;----------------------------------------------------------------------;

;------------  Téléchargement de données pour les graphes  ------------;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;----------------------------------------------------------------------;
begin
;------------------------  Lecture des champs  ------------------------;
NLEV  = numAsciiRow(dossier+"ZY"+sprinti("%0.11i",0))
NYPT  = numAsciiCol(dossier+"ZY"+sprinti("%0.11i",0))
z = asciiread(dossier+"ZY"+sprinti("%0.11i",0),(/NLEV,NYPT/),"double")
y = asciiread(dossier+"YZ"+sprinti("%0.11i",0),(/NLEV,NYPT/),"double")
filename = dossier+var+sprinti("%0.11i",step)
RES = asciiread(filename,(/NLEV,NYPT/),"double")

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
map@sfXArray             = y
map@vpWidthF             = 0.25
map@tiXAxisString       = "y(km)"
map@tiXAxisFontHeightF  = 0.005 
map@tmXBMinorPerMajor    = 1
map@tmXBMode             = "Manual"
map@tmXBLabelFontAspectF = 1.5
map@tmXBLabelFontHeightF = .005
map@trXMinF              = 0
map@tmXBTickStartF       = 0
map@tmXBTickSpacingF     = 1000
map@trXMaxF              = 6000
map@tmXBTickEndF         = 6000
map@tmXBMajorOutwardLengthF = 0.0
map@tmXBMinorOutwardLengthF = 0
map@tmXBMajorLengthF = 0.0025
map@tmXBMajorThicknessF = 0.5
map@tmXBMinorLengthF = 0


map@sfYArray             = 1-z
map@vpHeightF            = 0.155
;map@tiYAxisSide         = "Left"
map@tiYAxisString       = "z(km)"
map@tiYAxisFontHeightF  = 0.005
map@tmYLMinorPerMajor    = 1
map@tmYLMode             = "Manual"
map@tmYLLabelFontAspectF = 1.5
map@tmYLLabelFontHeightF = .005
map@trYMinF              = 0
map@tmYLTickStartF       = 0
map@tmYLTickSpacingF     = 0.1
map@trYMaxF              = 1
map@tmYLTickEndF         = 1
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
map@lbBoxLinesOn               = True
map@lbLabelAutoStride          = False
map@cnLabelBarEndStyle         = "ExcludeOuterBoxes"
map@pmLabelBarWidthF           = .02
map@pmLabelBarHeightF          = .165
map@lbLabelFontThicknessF      = 0.8
map@lbLabelFontHeightF         =.004

; Marque les ligne de séparation sans afficher la valeur
map@cnLinesOn           = True
map@cnLineLabelPerimOn  = False
map@cnLabelMasking      = False
map@cnLevelFlags        = False
map@cnLineThicknessF    = 0.5
map@cnLineLabelsOn      = False
map@gsnContourNegLineDashPattern = 1
;----------------------------------------------------------------------;

;------------------------  Choix des couleurs  ------------------------;
map@cnFillOn = False
wks  = gsn_open_wks("epsi",var+sprinti("%0.11i",step))

; Général de [-max ; +max]
  gsn_define_colormap(wks,"temp_19lev")
  map@cnLevelSelectionMode = "ManualLevels"
  map@cnMaxLevelValF  = doubletointeger(10*max((/abs(max(RES(:,:))),abs(min(RES(:,:)))/)))/10.
  ;map@cnMinLevelValF  = -map@cnMaxLevelValF
map@cnMinLevelValF  =doubletointeger(10*min((/abs(max(RES(:,:))),abs(min(RES(:,:)))/)))/10.
  map@cnLevelSpacingF = (map@cnMaxLevelValF-map@cnMinLevelValF)/15
  map@cnFillColors = (/ 2,3,4,5,6,7,8,9,0,14,15,16,17,18,19,20,21 /)


; HOT BUBLLE
;  gsn_define_colormap(wks,"temp_19lev")
;  map@cnLevelSelectionMode = "ExplicitLevels"
;  map@cnLevels     = (/ -0.003,-0.00275,-0.0025,-0.00225,-0.0020,-0.00175,-0.0015,-0.00125,-0.001,-0.00075,-0.0005,0.0005,0.00075,0.001,0.00125,0.0015,0.00175,0.002,0.00225,0.0025,0.00275,0.003 /)
;  map@cnFillColors = (/ 2,3,3,4,4,5,6,7,7,8,9,0,14,15,16,17,17,18,19,19,20,20,21 /)
;  map@cnExplicitLabelBarLabelsOn = True 
;  map@lbLabelStrings             = (/ "-3","-2.75","-2.5","-2.25","-2","-1.75","-1.5","-1.25","-1","-0.75","-0.5","0.5","0.75","1","1.25","1.5","1.75","2","2.25","2.5","2.75","3" /)

  ;gsn_define_colormap(wks,"MPL_Greys")
 ; map@cnLevelSelectionMode = "ExplicitLevels"
 ; map@cnLevels     = (/ -0.003,-0.00275,-0.0025,-0.00225,-0.0020,-0.00175,-0.0015,-0.00125,-0.001,-0.00075,-0.0005,0.0005,0.00075,0.001,0.00125,0.0015,0.00175,0.002,0.00225,0.0025,0.00275,0.003 /)
 ; map@cnFillColors = (/ 120,110,99,87,78,66,57,45,36,25,15,0,15,25,36,45,57,66,78,87,99,110,120 /)
 ; map@cnExplicitLabelBarLabelsOn = True 
 ; map@lbLabelStrings             = (/ "-3","-2.75","-2.5","-2.25","-2","-1.75","-1.5","-1.25","-1","-0.75","-0.5","0.5","0.75","1","1.25","1.5","1.75","2","2.25","2.5","2.75","3" /)


 ; gsn_define_colormap(wks,"MPL_Greys")
 ; map@cnLevelSelectionMode = "ExplicitLevels"
 ; map@cnLevels     = (/ -0.003,-0.0025,-0.0020,-0.0015,-0.001,-0.0005,0.,0.0006,0.001,0.0015,0.002,0.0025,0.003 /)
 ; map@cnFillColors = (/ 120,99,78,57,36,15,0,0,15,36,57,78,99,120 /)
 ; map@cnExplicitLabelBarLabelsOn = True 
 ; map@lbLabelStrings             = (/ "-3","-2.5","-2","-1.5","-1","-0.5","0.0","0.5","1","1.5","2","2.5","3" /)


;  gsn_define_colormap(wks,"MPL_Greys")
;  map@cnLevelSelectionMode = "ExplicitLevels"
;  map@cnLevels     = (/ 0.,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009 /)
;  map@cnFillColors = (/ 0,0,12,24,36,48,60,72,84,96,108,120 /)
;  map@cnExplicitLabelBarLabelsOn = True 
;  map@lbLabelStrings             = (/ "0","1","2","3","4","5","6","7","8","9","10" /)

;----------------------------------------------------------------------;
;-----------------------  Création des graphes  -----------------------;

plot = gsn_csm_contour(wks,RES(:,:),map)
; graphe final
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
