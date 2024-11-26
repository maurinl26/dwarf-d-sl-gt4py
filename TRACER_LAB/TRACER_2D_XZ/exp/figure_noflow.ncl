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
experience = "TEST"
dossier = "./"+experience+"/res/"
;dossier = "../res/"
;----------------------------------------------------------------------;

;------------  Téléchargement de données pour les graphes  ------------;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;----------------------------------------------------------------------;
;----------  Efface les messages d'erreur  ----------------------------;
;err = NhlGetErrorObjectId()
;setvalues err
;  "errLevel" : "Fatal"
;end setvalues
;----------------------------------------------------------------------;
;----------------------------------------------------------------------;
begin
;------------------------  Lecture des champs  ------------------------;
NLEV  = numAsciiRow(dossier+"Z"+sprinti("%0.11i",0))
NXPT  = numAsciiCol(dossier+"Z"+sprinti("%0.11i",0))
z = asciiread(dossier+"Z"+sprinti("%0.11i",0),(/NLEV,NXPT/),"double")
x = asciiread(dossier+"X"+sprinti("%0.11i",0),(/NLEV,NXPT/),"double")

print("")
print("DIMMENSIONS DU DOMAINE :")
print("NXPT = "+NXPT+" ; NLEV = "+NLEV)
print("")

filename = dossier+var+sprinti("%0.11i",step)
RES = asciiread(filename,(/NLEV,NXPT/),"double")
RES1tmp = asciiread(filename,(/NLEV,NXPT/),"double")

do i = 2, NXPT-3
  do j = 0, NLEV-1
     RES(j,i) = (6*RES1tmp(j,i)+4*RES1tmp(j,i+1)+4*RES1tmp(j,i-1)+RES1tmp(j,i+2)+RES1tmp(j,i-2))/16
;     RES(j,i) = (2*RES1tmp(j,i) + RES1tmp(j,i+1)+RES1tmp(j,i-1))/4
;     RES(j,i) = (RES1tmp(j,i+1)+RES1tmp(j,i-1))/2
;     RES(j,i) = RES1tmp(j,i)
 end do
end do

filename = dossier+"L"+sprinti("%0.11i",step)
TEMP = asciiread(filename,(/NLEV,NXPT/),"double")

;----------------------------------------------------------------------;

;----------------  Valeurs des paramètres des graphes  ----------------;
map  = True
mapT = True

; Enlève des infos inutiles
map@cnInfoLabelOn = False
map@gsnFrame = False
map@gsnDraw  = False

mapT@cnInfoLabelOn = False
mapT@gsnFrame = False
mapT@gsnDraw  = False

; Axes abscisse et ordonnée
map@cnFillMode = "AreaFill"
map@trGridType = "TriangularMesh"

map@sfXArray                = x
map@tiXAxisString           = "x(km)"
map@tiXAxisFontHeightF      = 0.0125
map@vpWidthF                = 0.6
map@tmXBMinorPerMajor       = 1
map@tmXBMode                = "Manual"
map@tmXBTickEndF            = max(x(:,:))
map@trXMaxF                 = max(x(:,:))
map@tmXBTickSpacingF        = 10
map@tmXBTickStartF          = min(x(:,:))
map@trXMinF                 = min(x(:,:))
map@tmXBLabelFontAspectF    = 1.5
map@tmXBLabelFontHeightF    = .01
map@tmXBMajorOutwardLengthF = 0.0
map@tmXBMinorOutwardLengthF = 0
map@tmXBMajorLengthF        = 0.0025
map@tmXBMajorThicknessF     = 0.5
map@tmXBMinorLengthF        = 0

map@sfYArray = z
map@tiYAxisSide             = "Left"
map@tiYAxisString           = "z(km)"
map@tiYAxisFontHeightF      = 0.0125
map@vpHeightF               = 0.6
map@tmYLMinorPerMajor       = 1
map@tmYLMode                = "Manual"
;map@tmYLTickEndF           = 20
map@trYMaxF                 = max(z(:,:))
map@tmYLTickEndF            = max(z(:,:))
map@tmYLTickSpacingF        = 1
map@trYMinF                 = min(z(:,:))
map@tmYLTickStartF          = 0
map@tmYLLabelFontAspectF    = 1.5
map@tmYLLabelFontHeightF    = .01
map@tmYLMajorOutwardLengthF = 0.0
map@tmYLMinorOutwardLengthF = 0
map@tmYLMajorLengthF        = 0.0025
map@tmYLMajorThicknessF     = 0.5
map@tmYLMinorLengthF        = 0
map@tmBorderThicknessF      = 0.8


mapT@sfXArray                = x
mapT@tiXAxisString           = "x(km)"
mapT@tiXAxisFontHeightF      = 0.0125
mapT@vpWidthF                = 0.6
mapT@tmXBMinorPerMajor       = 1
mapT@tmXBMode                = "Manual"
mapT@tmXBTickEndF            = max(x(:,:))
mapT@trXMaxF                 = max(x(:,:))
mapT@tmXBTickSpacingF        = 10
mapT@tmXBTickStartF          = min(x(:,:))
mapT@trXMinF                 = min(x(:,:))
mapT@tmXBLabelFontAspectF    = 1.5
mapT@tmXBLabelFontHeightF    = .01
mapT@tmXBMajorOutwardLengthF = 0.0
mapT@tmXBMinorOutwardLengthF = 0
mapT@tmXBMajorLengthF        = 0.0025
mapT@tmXBMajorThicknessF     = 0.5
mapT@tmXBMinorLengthF        = 0


mapT@sfYArray = z
;mapT@tiYAxisSide            = "Left"
mapT@tiYAxisString           = "z(km)"
mapT@tiYAxisFontHeightF      = 0.0125
mapT@vpHeightF               = 0.6
mapT@tmYLMinorPerMajor       = 1
mapT@tmYLMode                = "Manual"
mapT@trYMaxF                 = max(z(:,:))
mapT@tmYLTickEndF            = max(z(:,:))
mapT@tmYLTickSpacingF        = 1
mapT@trYMinF                 = min(z(:,:))
mapT@tmYLTickStartF          = 0
mapT@tmYLLabelFontAspectF    = 1.5
mapT@tmYLLabelFontHeightF    = .01
mapT@tmYLMajorOutwardLengthF = 0.0
mapT@tmYLMinorOutwardLengthF = 0
mapT@tmYLMajorLengthF        = 0.0025
mapT@tmYLMajorThicknessF     = 0.5
mapT@tmYLMinorLengthF        = 0
mapT@tmBorderThicknessF      = 0.8



; Légende
map@lbLabelBarOn               = True
map@lbOrientation              = "vertical" 
map@cnLineLabelBackgroundColor = "transparent"
map@lbBoxLinesOn               = False
map@lbLabelAutoStride          = False
map@cnLabelBarEndStyle         = "ExcludeOuterBoxes"
map@lbLabelFontHeightF         =.00925
;map@lbTitleOn                 = True
;map@lbTitleString             = "x10~S~-3"
;map@lbTitleFontHeightF        = 0.025
;map@pmLabelBarOrthogonalPosF  = .001 
;map@pmLabelBarParallelPosF    = 0.51
map@pmLabelBarWidthF           = .04
map@pmLabelBarHeightF          = .62


; Légende
mapT@lbLabelBarOn               = True
mapT@lbOrientation              = "vertical" 
mapT@cnLineLabelBackgroundColor = "transparent"
mapT@lbBoxLinesOn               = True
mapT@lbLabelAutoStride          = False
mapT@cnLabelBarEndStyle         = "ExcludeOuterBoxes"
mapT@pmLabelBarWidthF           = .02
mapT@pmLabelBarHeightF          = .165
mapT@lbLabelFontThicknessF      = 0.8
mapT@lbLabelFontHeightF         =.004


; Marque les ligne de séparation sans afficher la valeur
map@cnFillOn            = True
map@cnLinesOn           = True
map@cnLineLabelPerimOn  = False
map@cnLabelMasking      = False
map@cnLevelFlags        = True
map@cnLineThicknessF    = 0.5
map@cnLineLabelsOn      = False
map@gsnContourNegLineDashPattern  = 2
;map@gsnContourPosLineDashPattern = 0


mapT@cnLinesOn                 = True
mapT@cnLineLabelPerimOn        = False
mapT@cnLineLabelFont           = 21 
mapT@cnLineLabelFontHeightF    = 0.0025 
mapT@cnLineLabelFontThicknessF = 0.03 
mapT@cnLabelMasking            = True
mapT@cnLevelFlags              = True
mapT@cnLineThicknessF          = 0.25
mapT@cnLineLabelInterval       = 2  
mapT@cnLineLabelPlacementMode  = "computed"
mapT@cnLineLabelsOn            = True
;mapT@cnLineDashSegLenF        = 0.1
mapT@gsnContourNegLineDashPattern = 1
;mapT@gsnContourZeroLineThicknessF = 2.

;----------------------------------------------------------------------;
;------------------------  Choix des couleurs  ------------------------;
;----------------------------------------------------------------------;

;----------------------------------------------------------------------;
wks  = gsn_open_wks("epsi",var+sprinti("%0.11i",step))
;----------------------------------------------------------------------;

; Général de [-max ; +max]
gsn_define_colormap(wks,"temp_19lev")
map@cnLevelSelectionMode = "ManualLevels"
  map@cnMaxLevelValF  = doubletointeger(10*max((/abs(max(RES(:,:))),abs(min(RES(:,:)))/)))/10.
  map@cnMinLevelValF  = -map@cnMaxLevelValF
  map@cnLevelSpacingF = (map@cnMaxLevelValF-map@cnMinLevelValF)/15
  map@cnFillColors = (/ 2,3,4,5,6,7,8,9,0,14,15,16,17,18,19,20,21 /)

; SCHAR
;gsn_define_colormap(wks,"cmp_b2r")
;map@cnLevelSelectionMode = "ExplicitLevels"
;  map@cnLevels     = (/ -1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1 /)
;  map@cnFillOpacityF = 1.0
;  map@cnFillColors = (/ 2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,24,26,28,0,37,39,41,43,45,47,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65 /)
;  map@cnExplicitLabelBarLabelsOn = True 
;  map@lbLabelStrings             = (/"-1.00","-0.95","-0.90","-0.85","-0.80","-0.75","-0.70","-0.65","-0.60","-0.55","-0.50","-0.45","-0.40","-0.35","-0.30","-0.25","-0.20","-0.15","-0.10","-0.05","0.05","0.10","0.15","0.20","0.25","0.30","0.35","0.40","0.45","0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90","0.95","1.00"/)

;gsn_define_colormap(wks,"MPL_Greys")
;map@cnLevelSelectionMode = "ExplicitLevels"
;  map@cnLevels     = (/ -1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1 /)
;  map@cnFillOpacityF = 1.0
;  map@cnFillColors = (/ 128,122,116,110,104,98,92,86,80,74,68,62,56,50,44,38,32,26,20,14,0,0,14,20,26,32,38,44,50,56,62,68,74,80,86,92,98,104,110,116,122,124,126,128 /)
;  map@cnExplicitLabelBarLabelsOn = True 
;  map@lbLabelStrings  =(/"-1.0"," ","-0.9"," ","-0.8"," ","-0.7"," ","-0.6"," ","-0.5"," ","-0.4"," ","-0.3"," ","-0.2"," ","-0.1"," ","0"," ","0.1"," ","0.2"," ","0.3"," ","0.4"," ","0.5"," ","0.6"," ","0.7"," ","0.8"," ","0.9"," ","1.0"/)

;gsn_define_colormap(wks,"MPL_Greys")
;map@cnLevelSelectionMode = "ExplicitLevels"
;  map@cnLevels     = (/ -1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1 /)
;  map@cnFillOpacityF = 1.0
;  map@cnFillColors = (/ 128,116,104,92,80,68,56,44,32,20,0,0,20,32,44,56,68,80,92,104,116,124,128 /)
;  map@cnExplicitLabelBarLabelsOn = True 
;  map@lbLabelStrings  =(/"-1.0","-0.9","-0.8","-0.7","-0.6","-0.5","-0.4","-0.3","-0.2","-0.1","0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"/)


plot  = gsn_csm_contour(wks,RES(:,:),map)

;----------------------------------------------------------------------;
;--------------------  Paramétrisation du relief  ---------------------;
; Dessine le relief

xpoly = new(NXPT+3,"double")
ypoly = new(NXPT+3,"double")
do i = 0, NXPT-1
   xpoly(i) = x(NLEV-1,i) 
   ypoly(i) = z(NLEV-1,i)
end do
xpoly(NXPT)     = x(NLEV-1,NXPT-1)       
ypoly(NXPT)     = 0
xpoly(NXPT+1)   = x(NLEV-1,0)   
ypoly(NXPT+1)   = 0
xpoly(NXPT+2)   = x(NLEV-1,0)    
ypoly(NXPT+2)   = 0        
pres                  = True
;pres@gsFillColor     = "Grey"
pres@gsFillColor      = "Black"
pres@gsLineThicknessF = 2.0

; dessine le relief
dum0 = gsn_add_polygon(wks,plot,xpoly,ypoly,pres)

;----------------------------------------------------------------------;
;----------------------------  Iso-Theta  -----------------------------;

mapT@cnLevelSelectionMode = "ManualLevels"
mapT@cnMaxLevelValF  = doubletointeger(10*max(TEMP(:,:)))/10
mapT@cnMinLevelValF  = doubletointeger(10*min(TEMP(:,:)))/10
mapT@cnLevelSpacingF = (mapT@cnMaxLevelValF-mapT@cnMinLevelValF)/15

;plotT = gsn_csm_contour(wks,TEMP(:,:),mapT)
;overlay(plot,plotT)

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
