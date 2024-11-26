;####################  SCRIPT NCL POUR TRACER DES CHAMPS CALCULÉS  #####################;
;#                                                                                     #;
;# auteurs : F.Voitus                              #;
;# sujet   : Permet de visualiser les résultats du dossier ${experience} via           #;
;#           la commande suivante :                                                    #;
;#                            ncl 'var="X"' 'step=n' figure.ncl                        #;
;#           où ${X} une des variables {U,W,T}                                         #;
;#                                                                                     #;
;#######################################################################################;

;---------------------  Definition des varialbes  ---------------------;
;experience = "HE_baroc_64Niv"
experience = "EE_baroc_64Niv"
;experience = "QE_baroc_64Niv"
;experience = "TEST"
dossier = "./"+experience+"/res/Z/"
;dossier = "../res/"
;----------------------------------------------------------------------;

x_max = 180
x_min  = -60
x_unite = 0.5
y_max = 3.5
y_min  =-3.5
y_unite = 0.5


;------------  Téléchargement de données pour les graphes  ------------;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;--------------------------------------------------------------------------------------;

;----------  Efface les messages d'erreur  -----------;
;err = NhlGetErrorObjectId()
;setvalues err
;  "errLevel" : "Fatal"
;end setvalues
;-------------------------------------------------------------;

begin

;------------------------  Lecture des champs  ------------------------;
NYPT  = numAsciiRow(dossier+"XY"+sprinti("%0.11i",0))
NXPT  = numAsciiCol(dossier+"XY"+sprinti("%0.11i",0))
y = asciiread(dossier+"YX"+sprinti("%0.11i",0),(/NYPT,NXPT/),"double")
x = asciiread(dossier+"XY"+sprinti("%0.11i",0),(/NYPT,NXPT/),"double")

filename = dossier+"TPXY"+sprinti("%0.11i",step)
TEMP = asciiread(filename,(/NYPT,NXPT/),"double")

filename = dossier+"PPXY"+sprinti("%0.11i",step)
PRES = asciiread(filename,(/NYPT,NXPT/),"double")
;PRES=PRES/100

print("")
print("DIMMENSIONS DU DOMAINE :")
print("NXPT = "+NXPT+" ; NYPT = "+NYPT)
print("")
;----------------------------------------------------------------------;
;----------------  Valeurs des paramètres des graphes  ----------------;
resT = True
resP = True
 
; Enlève des infos inutiles
resP@cnInfoLabelOn = False
resP@gsnFrame = False
resP@gsnDraw  = False

resT@cnInfoLabelOn = False
resT@gsnFrame = False
resT@gsnDraw  = False

; Axes abscisse et ordonnée
resP@cnFillMode = "AreaFill"
resP@trGridType = "TriangularMesh"
resP@sfXArray             = x/100
resP@vpWidthF             = 0.2
resP@tiXAxisString       = "x(100km)"
resP@tiXAxisFontHeightF  = 0.003
resP@tmXBMinorPerMajor    = 1
resP@tmXBMode             = "Manual"
resP@tmXBLabelFontAspectF = 1.5
resP@tmXBLabelFontHeightF = .003
resP@trXMinF              = 25
resP@tmXBTickStartF       = 25
resP@tmXBTickSpacingF     = 20
resP@trXMaxF              = 150
resP@tmXBTickEndF         = 150
resP@tmXBMajorOutwardLengthF = 0.0
resP@tmXBMinorOutwardLengthF = 0
resP@tmXBMajorLengthF = 0.001
resP@tmXBMajorThicknessF = 0.25
resP@tmXBMinorLengthF = 0

resP@sfYArray             = y/100
resP@vpHeightF            = 0.1
;resP@tiYAxisSide         = "Left"
resP@tiYAxisString       = "y(100km)"
resP@tiYAxisFontHeightF  = 0.003
resP@tmYLMinorPerMajor    = 1
resP@tmYLMode             = "Manual"
resP@tmYLLabelFontAspectF = 1.5
resP@tmYLLabelFontHeightF = .003
resP@trYMinF              = 0
resP@tmYLTickStartF       = 0
resP@tmYLTickSpacingF     = 10
resP@trYMaxF              = 60
resP@tmYLTickEndF         = 61
resP@tmYLMajorOutwardLengthF = 0.0
resP@tmYLMinorOutwardLengthF = 0.0
resP@tmYLMajorLengthF = 0.001
resP@tmYLMajorThicknessF = 0.25
resP@tmYLMinorLengthF = 0
resP@tmBorderThicknessF = 0.8


resT@cnFillMode = "AreaFill"
resT@trGridType = "TriangularMesh"
resT@sfXArray             = x/100
resT@vpWidthF             = 0.2
resT@tiXAxisString       = "x(100km)"
resT@tiXAxisFontHeightF  = 0.003
resT@tmXBMinorPerMajor    = 1
resT@tmXBMode             = "Manual"
resT@tmXBLabelFontAspectF = 1.5
resT@tmXBLabelFontHeightF = .003
resT@trXMinF              = 25
resT@tmXBTickStartF       = 25
resT@tmXBTickSpacingF     = 20
resT@trXMaxF              = 150
resT@tmXBTickEndF         = 150
resT@tmXBMajorOutwardLengthF = 0.0
resT@tmXBMinorOutwardLengthF = 0
resT@tmXBMajorLengthF = 0.001
resT@tmXBMajorThicknessF = 0.25
resT@tmXBMinorLengthF = 0

resT@sfYArray             = y/100
resT@vpHeightF            = 0.1
resT@tiYAxisString       = "y(100km)"
resT@tiYAxisFontHeightF  = 0.003
resT@tmYLMinorPerMajor    = 1
resT@tmYLMode             = "Manual"
resT@tmYLLabelFontAspectF = 1.5
resT@tmYLLabelFontHeightF = .003
resT@trYMinF              = 0
resT@tmYLTickStartF       = 0
resT@tmYLTickSpacingF     = 10
resT@trYMaxF              = 60
resT@tmYLTickEndF         = 61
resT@tmYLMajorOutwardLengthF = 0.0
resT@tmYLMinorOutwardLengthF = 0.0
resT@tmYLMajorLengthF = 0.001
resT@tmYLMajorThicknessF = 0.25
resT@tmYLMinorLengthF = 0
resT@tmBorderThicknessF = 0.8

; Légende
resP@lbLabelBarOn               = True
resP@lbOrientation              = "vertical" 
resP@cnLineLabelBackgroundColor = "transparent"
resP@lbBoxLinesOn               = True
resP@lbLabelAutoStride          = False
resP@cnLabelBarEndStyle         = "ExcludeOuterBoxes"
resP@pmLabelBarWidthF           = .015
resP@pmLabelBarHeightF          = .105
resP@lbLabelFontThicknessF      = 0.8
resP@lbLabelFontHeightF         =.003

; Légende
resT@lbLabelBarOn               = False
resT@lbOrientation              = "vertical" 
resT@cnLineLabelBackgroundColor = "transparent"
resT@lbBoxLinesOn               = True
resT@lbLabelAutoStride          = False
resT@cnLabelBarEndStyle         = "ExcludeOuterBoxes"
resT@pmLabelBarWidthF           = .02
resT@pmLabelBarHeightF          = .165
resT@lbLabelFontThicknessF      = 0.8
resT@lbLabelFontHeightF         =.004

; Marque les ligne de séparation sans afficher la valeur
resP@cnLinesOn           = False
;resP@cnLineLabelPerimOn  = False
;resP@cnLineLabelFont = 21 
;resP@cnLineLabelFontHeightF = 0.0025 
;resP@cnLineLabelFontThicknessF = 0.03 
;resP@cnLabelMasking      = True
;resP@cnLevelFlags        = True
;resP@cnLineThicknessF    = 0.25
;resP@cnLineLabelInterval   = 2  
;resP@cnLineLabelPlacementMode = "computed"
;resP@cnLineLabelsOn      = False



;Marque les ligne de séparation sans afficher la valeur
resT@cnLinesOn           = True
resT@cnLineLabelPerimOn  = False
resT@cnLineLabelFont = 21 
resT@cnLineLabelFontHeightF = 0.0025 
resT@cnLineLabelFontThicknessF = 0.03 
resT@cnLabelMasking      = True
resT@cnLevelFlags        = True
resT@cnLineThicknessF    = 0.25
resT@cnLineLabelInterval   = 2  
resT@cnLineLabelPlacementMode = "computed"
resT@cnLineLabelsOn      = True
;resT@cnLineDashSegLenF     = 0.1
resT@gsnContourNegLineDashPattern = 1
;resT@gsnContourZeroLineThicknessF = 2.




;------------------------  Choix des couleurs  ------------------------;
resP@cnFillOn = True

wks  = gsn_open_wks("epsi",experience+"_10days")

resP@cnLevelSelectionMode = "ManualLevels"
resT@cnLevelSelectionMode = "ManualLevels"

resP@cnMaxLevelValF  = doubletointeger(10*960)/10
resP@cnMinLevelValF  = doubletointeger(10*930)/10
resP@cnLevelSpacingF = (resP@cnMaxLevelValF-resP@cnMinLevelValF)/15

resT@cnMaxLevelValF  = doubletointeger(10*300)/10
resT@cnMinLevelValF  = doubletointeger(10*268)/10
resT@cnLevelSpacingF = (resT@cnMaxLevelValF-resT@cnMinLevelValF)/15

   

gsn_define_colormap(wks,"MPL_Greys")
resP@cnLevelSelectionMode = "ExplicitLevels"
resP@cnLevels     = (/ 930,932,934,936,938,940,942,944,946,948,950,952,954,956,958,960/)
 
;----------------------------------------------------------------------;
;-----------------------  Création des graphes  -----------------------;

plot_resP = gsn_csm_contour(wks,PRES(:,:),resP)
plot_resT = gsn_csm_contour(wks,TEMP(:,:),resT)
overlay(plot_resP,plot_resT)

; graphe final
draw(plot_resP)
frame(wks)
delete(wks)
;----------------------------------------------------------------------;

;----------------------  Visionnage des graphes  ----------------------;
system("epstopdf "+experience+"_10days"+".epsi")
system("evince "+experience+"_10days"+".pdf")
system("rm -f *.epsi *.pdf")
system("ls -t "+dossier+var+"*")
;----------------------------------------------------------------------;

end
