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
;experience = "HE_baroc"
experience = "EE_baroc15"
;experience = "QE_baroc"
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
PRES=PRES/100

print("")
print("DIMMENSIONS DU DOMAINE :")
print("NXPT = "+NXPT+" ; NYPT = "+NYPT)
print("")
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
map@sfXArray             = x/100
map@vpWidthF             = 0.25
map@tiXAxisString       = "~F14~x~F13~(100km)"
map@tiXAxisFontHeightF  = 0.004
map@tmXBMinorPerMajor    = 1
map@tmXBMode             = "Manual"
map@tmXBLabelFontAspectF = 1.5
map@tmXBLabelFontHeightF = .004
map@trXMinF              = 0
map@tmXBTickStartF       = 0
map@tmXBTickSpacingF     = 50
map@trXMaxF              = 200
map@tmXBTickEndF         = 200
map@tmXBMajorOutwardLengthF = 0.0
map@tmXBMinorOutwardLengthF = 0
map@tmXBMajorLengthF = 0.0025
map@tmXBMajorThicknessF = 0.5
map@tmXBMinorLengthF = 0


map@sfYArray             = y/100
map@vpHeightF            = 0.155
;map@tiYAxisSide         = "Left"
map@tiYAxisString       = "~F14~y~F13~(100km)"
map@tiYAxisFontHeightF  = 0.004
map@tmYLMinorPerMajor    = 1
map@tmYLMode             = "Manual"
map@tmYLLabelFontAspectF = 1.5
map@tmYLLabelFontHeightF = .004
map@trYMinF              = 0
map@tmYLTickStartF       = 0
map@tmYLTickSpacingF     = 10
map@trYMaxF              = 60
map@tmYLTickEndF         = 61
map@tmYLMajorOutwardLengthF = 0.0
map@tmYLMinorOutwardLengthF = 0.0
map@tmYLMajorLengthF = 0.0025
map@tmYLMajorThicknessF = 0.5
map@tmYLMinorLengthF = 0

map@tmBorderThicknessF = 0.8

; Légende
map@lbLabelBarOn               = False
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



;----------------------------------------------------------------------;

;------------------------  Choix des couleurs  ------------------------;
map@cnFillOn = False
wks  = gsn_open_wks("epsi","baro"+sprinti("%0.11i",step))
gsn_define_colormap(wks,"BlueWhiteOrangeRed")

map@cnLinesOn           = True
map@cnLineLabelPerimOn  = False
map@cnLineLabelFont = 21 
map@cnLineLabelFontHeightF = 0.0025 
map@cnLineLabelFontThicknessF = 0.03 
map@cnLabelMasking      = True
map@cnLevelFlags        = True
map@cnLineThicknessF    = 0.4
map@cnLineLabelInterval   = 2  
;map@gsnContourPosLineDashPattern = 2


; Général de [-max ; +max]
map@cnLevelSelectionMode = "ManualLevels"
map@cnMaxLevelValF  = doubletointeger(10*302)/10
map@cnMinLevelValF  =doubletointeger(10*266)/10
map@cnLevelSpacingF = (map@cnMaxLevelValF-map@cnMinLevelValF)/18
map@cnLineLabelPlacementMode = "computed"
map@cnLineColors = "Grey"
 

 
;----------------------------------------------------------------------;
;-----------------------  Création des graphes  -----------------------;

plot_temp = gsn_csm_contour(wks,TEMP(:,:),map)



map@cnMaxLevelValF  = doubletointeger(10*1000)/10
map@cnMinLevelValF  =doubletointeger(10*930)/10
map@cnLevelSpacingF = (map@cnMaxLevelValF-map@cnMinLevelValF)/50
map@cnLineLabelPlacementMode = "constant"
map@cnLineColors = 1
 map@gsnContourPosLineDashPattern = 1

plot_pres = gsn_csm_contour(wks,PRES(:,:),map)

overlay(plot_temp,plot_pres)

; graphe final
draw(plot_temp)
frame(wks)
delete(wks)
;----------------------------------------------------------------------;

;----------------------  Visionnage des graphes  ----------------------;
system("epstopdf "+"baro"+sprinti("%0.11i",step)+".epsi")
system("evince "+"baro"+sprinti("%0.11i",step)+".pdf")
system("rm -f *.epsi *.pdf")
;system("ls -t "+dossier+var+"*")
;----------------------------------------------------------------------;

end
