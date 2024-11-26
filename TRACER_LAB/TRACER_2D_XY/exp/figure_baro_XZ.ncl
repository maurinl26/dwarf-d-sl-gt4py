;####################  SCRIPT NCL POUR TRACER DES CHAMPS CALCULÉS  #####################;
;#                                                                                     #;
;# auteurs : L.Auger, F.Voitus, C.Colavolpe, B.Menetrier                               #;
;# sujet   : Permet de visualiser les résultats du dossier ${experience} via           #;
;#           la commande suivante :                                                    #;
;#                            ncl 'var="X"' 'step=n' figure_baro_XZ.ncl                        #;
;#           où ${X} une des variables {U,W,T}                                         #;
;#                                                                                     #;
;#######################################################################################;

;---------------------  Definition des varialbes  ---------------------;
;experience = "EE_bplan_XZsect"
experience = "QE_bplan_XZsect"
dossier = "./"+experience+"/res/Y/"
;dossier = "../res/"
;----------------------------------------------------------------------;

;------------  Téléchargement de données pour les graphes  ------------;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;----------------------------------------------------------------------;
begin
;------------------------  Lecture des champs  ------------------------;
NLEV  = numAsciiRow(dossier+"ZX"+sprinti("%0.11i",0))
NXPT  = numAsciiCol(dossier+"ZX"+sprinti("%0.11i",0))
z = asciiread(dossier+"ZX"+sprinti("%0.11i",0),(/NLEV,NXPT/),"double")
x = asciiread(dossier+"XZ"+sprinti("%0.11i",0),(/NLEV,NXPT/),"double")
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
map@sfXArray             = x/100
map@vpWidthF             = 0.25
map@tiXAxisString       = " "
map@tiXAxisFontHeightF  = 0.003
map@tmXBMinorPerMajor    = 1

map@tmXBMode             = "Manual"
map@tmXBLabelFontAspectF = 1.5
map@tmXBLabelFontHeightF = .003
map@trXMinF              = 15
map@tmXBTickStartF       = 15
map@tmXBTickSpacingF     = 20
map@trXMaxF              = 154
map@tmXBTickEndF         = 154
map@tmXBMajorOutwardLengthF = 0.0
map@tmXBMinorOutwardLengthF = 0
map@tmXBMajorLengthF = 0.001
map@tmXBMajorThicknessF = 0.25
map@tmXBMinorLengthF = 0


map@sfYArray             = 1-z
map@vpHeightF            = 0.08
;map@tiYAxisSide         = "Left"
map@tiYAxisString       = "~F33~s"
map@tiYAxisFontHeightF  = 0.003
map@tmYLMinorPerMajor    = 1
map@tmYLMode             = "Manual"
map@tmYLLabelFontAspectF = 1.5
map@tmYLLabelFontHeightF = .003
map@trYMinF              = 0
map@tmYLTickStartF       = 0
map@tmYLTickSpacingF     = 0.1
map@trYMaxF              = 1
map@tmYLTickEndF         = 1
map@tmYLMajorOutwardLengthF = 0.0
map@tmYLMinorOutwardLengthF = 0.0
map@tmYLMajorLengthF = 0.002
map@tmYLMajorThicknessF = 0.25
map@tmYLMinorLengthF = 0.001
map@tmYLMinorThicknessF = 0.25
map@tmBorderThicknessF = 0.8

map@tmYLAutoPrecision = False
map@tmYLPrecision = 1

;map@tmYLMode   = "Explicit"                          
;map@tmYLValues = (/0.0078,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99/)
;map@tmYLLabelDeltaF=-0.9
;map@tmYLLabels = (/"0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1."/)
;map@tmYLMinorValues = (/"-80","-70","-50","-40","-20","-10","10","20","40","50","70","80"/)


; Légende
map@lbLabelBarOn               = True
map@lbOrientation              = "horizontal" 
map@cnLineLabelBackgroundColor = "transparent"
map@lbBoxLinesOn               = True
map@lbLabelAutoStride          = False
map@cnLabelBarEndStyle         = "ExcludeOuterBoxes"
map@pmLabelBarWidthF           = .2
map@pmLabelBarHeightF          = .015
map@lbLabelFontThicknessF      = 0.5
map@lbLabelFontHeightF         =.003

; Marque les ligne de séparation sans afficher la valeur
map@cnLinesOn           = True
map@cnLineLabelPerimOn  = False
map@cnLabelMasking      = False
map@cnLevelFlags        = False
map@cnLineThicknessF    = 0.2
map@cnLineLabelsOn      = False
map@gsnContourNegLineDashPattern = 2
;----------------------------------------------------------------------;

;------------------------  Choix des couleurs  ------------------------;
map@cnFillOn = True
wks  = gsn_open_wks("epsi",experience+"_"+sprinti("%0.11i",step))

; Général de [-max ; +max]
;  gsn_define_colormap(wks,"temp_19lev")
;  map@cnLevelSelectionMode = "ManualLevels"
;  map@cnMaxLevelValF  = doubletointeger(10*max((/abs(max(RES(:,:))),abs(min(RES(:,:)))/)))/10.
;  map@cnMinLevelValF  = -map@cnMaxLevelValF
;  map@cnLevelSpacingF = (map@cnMaxLevelValF-map@cnMinLevelValF)/15
;  map@cnFillColors = (/ 2,3,4,5,6,7,8,9,0,14,15,16,17,18,19,20,21 /)



 gsn_define_colormap(wks,"MPL_Greys")
  map@cnLevelSelectionMode = "ExplicitLevels"
  map@cnLevels     = (/ -17.,-15.,-13.,-11.,-9,-7.,-5.,-3.,-1.,1,3.,5.,7.,9.,11.,13.,15.,17./)
  map@cnFillColors = (/ 120,106,92,78,64,50,36,22,10,0,10,22,36,50,64,78,92,106,120/)
  map@cnExplicitLabelBarLabelsOn = True 
  map@lbLabelStrings             = (/ "-17","-15","-13","-11","-9","-7","-5","-3","-1","1","3","5","7","9","11","13","15" ,"17"/)


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
system("epstopdf "+experience+"_"+sprinti("%0.11i",step)+".epsi")
system("evince "+experience+"_"+sprinti("%0.11i",step)+".pdf")
system("rm -f *.epsi *.pdf")
system("ls -t "+dossier+var+"*")
;----------------------------------------------------------------------;

end
