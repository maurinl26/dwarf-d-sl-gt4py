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
experience = "EE_baroc"
;experience = "QE_baroc"
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
filename = dossier+var+sprinti("%0.11i",step)
RES = asciiread(filename,(/NYPT,NXPT/),"double")
print("")
print("DIMMENSIONS DU DOMAINE :")
print("NXPT = "+NXPT+" ; NYPT = "+NYPT)
print("")
;----------------------------------------------------------------------;

RES = RES/100

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
map@vpWidthF             = 0.2
map@tiXAxisString       = "x(100km)"
map@tiXAxisFontHeightF  = 0.003
map@tmXBMinorPerMajor    = 1

map@tmXBMode             = "Manual"
map@tmXBLabelFontAspectF = 1.5
map@tmXBLabelFontHeightF = .003
map@trXMinF              = 25
map@tmXBTickStartF       = 25
map@tmXBTickSpacingF     = 20
map@trXMaxF              = 150
map@tmXBTickEndF         = 150
map@tmXBMajorOutwardLengthF = 0.0
map@tmXBMinorOutwardLengthF = 0
map@tmXBMajorLengthF = 0.001
map@tmXBMajorThicknessF = 0.25
map@tmXBMinorLengthF = 0


map@sfYArray             = y/100
map@vpHeightF            = 0.1
;map@tiYAxisSide         = "Left"
map@tiYAxisString       = "y(100km)"
map@tiYAxisFontHeightF  = 0.003
map@tmYLMinorPerMajor    = 1
map@tmYLMode             = "Manual"
map@tmYLLabelFontAspectF = 1.5
map@tmYLLabelFontHeightF = .003
map@trYMinF              = 0
map@tmYLTickStartF       = 0
map@tmYLTickSpacingF     = 10
map@trYMaxF              = 60
map@tmYLTickEndF         = 61
map@tmYLMajorOutwardLengthF = 0.0
map@tmYLMinorOutwardLengthF = 0.0
map@tmYLMajorLengthF = 0.001
map@tmYLMajorThicknessF = 0.25
map@tmYLMinorLengthF = 0

map@tmBorderThicknessF = 0.8

; Légende
map@lbLabelBarOn               = True
map@lbOrientation              = "vertical" 
map@cnLineLabelBackgroundColor = "transparent"
map@lbBoxLinesOn               = True
map@lbLabelAutoStride          = False
map@cnLabelBarEndStyle         = "ExcludeOuterBoxes"
map@pmLabelBarWidthF           = .015
map@pmLabelBarHeightF          = .105
map@lbLabelFontThicknessF      = 0.8
map@lbLabelFontHeightF         =.003

; Marque les ligne de séparation sans afficher la valeur
map@cnLinesOn           = False
map@cnLineLabelPerimOn  = False
map@cnLineLabelFont = 21 
map@cnLineLabelFontHeightF = 0.0025 
map@cnLineLabelFontThicknessF = 0.03 
map@cnLabelMasking      = True
map@cnLevelFlags        = True
map@cnLineThicknessF    = 0.25
map@cnLineLabelInterval   = 2  
map@cnLineLabelPlacementMode = "computed"
map@cnLineLabelsOn      = False
;map@cnLineDashSegLenF     = 0.1
map@gsnContourNegLineDashPattern = 1
;map@gsnContourZeroLineThicknessF = 2.




;title figure
;map@tiMainString               = "Randomized: Default"
;----------------------------------------------------------------------;

;------------------------  Choix des couleurs  ------------------------;
map@cnFillOn = True
wks  = gsn_open_wks("epsi",var+sprinti("%0.11i",step))


map@cnLevelSelectionMode = "ManualLevels"

if (var .eq. "TPXY") then
  map@cnMaxLevelValF  = doubletointeger(10*300)/10
  map@cnMinLevelValF  =doubletointeger(10*268)/10
  map@cnLevelSpacingF = (map@cnMaxLevelValF-map@cnMinLevelValF)/15
else if (var .eq. "PPXY") then
  map@cnMaxLevelValF  = doubletointeger(10*960)/10
  map@cnMinLevelValF  = doubletointeger(10*930)/10
  map@cnLevelSpacingF = (map@cnMaxLevelValF-map@cnMinLevelValF)/15
else if (var .eq. "PSXY") then
  map@cnMaxLevelValF  = doubletointeger(10*1025)/10
  map@cnMinLevelValF  =doubletointeger(10*975)/10
  map@cnLevelSpacingF = (map@cnMaxLevelValF-map@cnMinLevelValF)/25
else
  map@cnMaxLevelValF  = doubletointeger(10*max((/abs(max(RES(:,:))),abs(min(RES(:,:)))/)))/10.
  map@cnMinLevelValF   =-doubletointeger(10*min((/abs(max(RES(:,:))),abs(min(RES(:,:)))/)))/10.
  map@cnLevelSpacingF = (map@cnMaxLevelValF-map@cnMinLevelValF)/15
end if
end if
end if

gsn_define_colormap(wks,"temp_19lev")
;gsn_define_colormap(wks,"MPL_Greys")
;map@cnLevelSelectionMode = "ExplicitLevels"
;map@cnLevels     = (/ 266,268,270,272,274,276,278,280,282,284,286,288,290,292,294,296,298,300,302/)
;map@cnMonoLineColor = False 
;map@cnLineColors = (/ 80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80/)
 
;gsn_define_colormap(wks,"temp_19lev")
gsn_define_colormap(wks,"MPL_Greys")
map@cnLevelSelectionMode = "ExplicitLevels"
map@cnLevels     = (/ 930,932,934,936,938,940,942,944,946,948,950,952,954,956,958,960/)
;map@cnMonoLineColor = False 
;map@cnLineColors = (/ 80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80/)
 
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
