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
dossier = "./"+experience+"/res/XZ/"
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
RES1tmp = asciiread(filename,(/NLEV,NXPT/),"double")

do i = 2, NXPT-3
  do j = 0, NLEV-1
;    RES(j,i) = (6*RES1tmp(j,i) + 4*RES1tmp(j,i+1)+4*RES1tmp(j,i-1)+RES1tmp(j,i+2)+RES1tmp(j,i-2))/16
;   RES(j,i) = (2*RES1tmp(j,i) + RES1tmp(j,i+1)+RES1tmp(j,i-1))/4
;   RES(j,i) = (RES1tmp(j,i+1)+RES1tmp(j,i-1))/2
   RES(j,i) = RES1tmp(j,i)
 end do
end do


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
map@tiXAxisString       = "x(km)"
map@tiXAxisFontHeightF  = 0.0125 
map@tmXBMinorPerMajor    = 1
map@tmXBMode             = "Manual"
map@tmXBLabelFontAspectF = 1.5
map@tmXBLabelFontHeightF = .005
map@trXMinF              = -20
map@tmXBTickStartF       = -20
map@tmXBTickSpacingF     = 5
map@trXMaxF              = 20
map@tmXBTickEndF         = 20
map@tmXBLabelFontAspectF    = 1.5
map@tmXBLabelFontHeightF    = .01
map@tmXBMajorOutwardLengthF = 0.0
map@tmXBMinorOutwardLengthF = 0
map@tmXBMajorLengthF        = 0.0025
map@tmXBMajorThicknessF     = 0.5
map@tmXBMinorLengthF        = 0


map@sfYArray               = z
map@tiYAxisFontHeightF = 0.0125
map@vpHeightF          = 0.6
;map@tiYAxisSide         = "Left"
map@tiYAxisString        = "z(km)"
map@tmYLMinorPerMajor  = 1
map@tmYLMode                 = "Manual"
map@tmYLLabelFontAspectF = 1.5
map@tmYLLabelFontHeightF = .005
map@trYMinF                      = 0
map@tmYLTickStartF          = 1
map@tmYLTickSpacingF     = 1
map@trYMaxF                 =  12
map@tmYLTickEndF            =  12
map@tmYLLabelFontAspectF    = 1.5
map@tmYLLabelFontHeightF    = .01
map@tmYLMajorOutwardLengthF = 0.0
map@tmYLMinorOutwardLengthF = 0
map@tmYLMajorLengthF        = 0.0025
map@tmYLMajorThicknessF     = 0.5
map@tmYLMinorLengthF        = 0
map@tmBorderThicknessF      = 0.8

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
map@cnFillOn = True
wks  = gsn_open_wks("epsi",var+sprinti("%0.11i",step))

; Général de [-max ; +max]
  gsn_define_colormap(wks,"temp_19lev")
  map@cnLevelSelectionMode = "ManualLevels"
  map@cnMaxLevelValF  = doubletointeger(10*max((/abs(max(RES(:,:))),abs(min(RES(:,:)))/)))/10.
  map@cnMinLevelValF  = -map@cnMaxLevelValF
  map@cnLevelSpacingF = (map@cnMaxLevelValF-map@cnMinLevelValF)/15
  map@cnFillColors = (/ 2,3,4,5,6,7,8,9,0,14,15,16,17,18,19,20,21 /)


; HOT BUBLLE
 ; gsn_define_colormap(wks,"temp_19lev")
 ; map@cnLevelSelectionMode = "ExplicitLevels"
 ; map@cnLevels     = (/ -0.003,-0.00275,-0.0025,-0.00225,-0.0020,-0.00175,-0.0015,-0.00125,-0.001,-0.00075,-0.0005,0.0005,0.00075,0.001,0.00125,0.0015,0.00175,0.002,0.00225,0.0025,0.00275,0.003 /)
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


; SCHAR
;gsn_define_colormap(wks,"cmp_b2r")
;map@cnLevelSelectionMode = "ExplicitLevels"
;  map@cnLevels     = (/ -1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1 /)
;  map@cnFillOpacityF = 1.0
;  map@cnFillColors = (/ 2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,24,26,28,0,37,39,41,43,45,47,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65 /)
;  map@cnExplicitLabelBarLabelsOn = True 
;  map@lbLabelStrings             = (/"-1.00","-0.95","-0.90","-0.85","-0.80","-0.75","-0.70","-0.65","-0.60","-0.55","-0.50","-0.45","-0.40","-0.35","-0.30","-0.25","-0.20","-0.15","-0.10","-0.05","0.05","0.10","0.15","0.20","0.25","0.30","0.35","0.40","0.45","0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90","0.95","1.00"/)



;  gsn_define_colormap(wks,"MPL_Greys")
;  map@cnLevelSelectionMode = "ExplicitLevels"
;  map@cnLevels     = (/ 0.,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009 /)
;  map@cnFillColors = (/ 0,0,12,24,36,48,60,72,84,96,108,120 /)
;  map@cnExplicitLabelBarLabelsOn = True 
;  map@lbLabelStrings             = (/ "0","1","2","3","4","5","6","7","8","9","10" /)

;----------------------------------------------------------------------;
;-----------------------  Création des graphes  -----------------------;

plot = gsn_csm_contour(wks,RES(:,:),map)

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
xpoly(NXPT+2)   = x(NLEV-1,NXPT/2+1)    
ypoly(NXPT+2)   = 0        
pres                  = True
;pres@gsFillColor     = "lightsalmon4"
pres@gsFillColor      = "Black"
pres@gsLineThicknessF = 2.0

; dessine le relief
dum0 = gsn_add_polygon(wks,plot,xpoly,ypoly,pres)

;----------------------------------------------------------------------;
;----------------------------  Iso-Theta  -----------------------------;





;----------------------------------------------------------------------;
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
