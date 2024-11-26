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
experience = "ARK_06_1"
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
map@vpWidthF             = 0.6
;map@tiXAxisString       = "Distance relative (en km)"
;map@tiXAxisFontHeightF  = 0.018
map@tmXBMinorPerMajor    = 1
map@tmXBMode             = "Manual"
map@tmXBLabelFontAspectF = 1.5
map@tmXBLabelFontHeightF = .015
map@trXMinF              = -25
map@tmXBTickStartF       = -25
map@tmXBTickSpacingF     = 50
map@trXMaxF              = 25
map@tmXBTickEndF         = 25

map@sfYArray             = z
map@vpHeightF            = 0.6
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
map@trYMaxF              = 12
map@tmYLTickEndF         = 12

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
map@pmLabelBarWidthF           = 0.04
map@pmLabelBarHeightF          = 0.6

; Marque les ligne de séparation sans afficher la valeur
map@cnLinesOn           = False
map@cnLineLabelPerimOn  = False
map@cnLabelMasking      = False
map@cnLevelFlags        = True
map@cnLineThicknessF    = 1
map@cnLineLabelsOn      = False
map@gsnContourNegLineDashPattern = 1
;----------------------------------------------------------------------;

;------------------------  Choix des couleurs  ------------------------;
map@cnFillOn = True
wks  = gsn_open_wks("epsi",var+sprinti("%0.11i",step))
gsn_define_colormap(wks,"temp_19lev")
;gsn_define_colormap(wks,"MPL_ocean")
;gsn_define_colormap(wks,"MPL_Greys")
;gsn_define_colormap(wks,"cmp_b2r")

; Général de [-max ; +max]
map@cnLevelSelectionMode = "ManualLevels"
  map@cnMaxLevelValF  = doubletointeger(10*max((/abs(max(RES(:,:))),abs(min(RES(:,:)))/)))/10.
  map@cnMinLevelValF  = -map@cnMaxLevelValF
  map@cnLevelSpacingF = (map@cnMaxLevelValF-map@cnMinLevelValF)/15
  map@cnFillColors = (/ 2,3,4,5,6,7,8,9,0,14,15,16,17,18,19,20,21 /)

;map@cnLevelSelectionMode = "ExplicitLevels"
  ;map@cnLevels     = (/ -1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1 /)
   ;map@cnLevels     = (/ -3,-2.85,-2.7,-2.55,-2.4,-2.25,-2.1,-1.95,-1.8,-1.65,-1.5,-1.35,-1.2,-1.05,-0.9,-0.75,-0.6,-0.45,-0.3,-0.15,0.15,0.3,0.45,0.6,0.75,0.9,1.05,1.2,1.35,1.5,1.65,1.8,1.95,2.1,2.25,2.4,2.55,2.7,2.85,3 /)
  ;map@cnExplicitLabelBarLabelsOn = False 
  ;map@cnFillOpacityF = 1.0
  ;map@cnFillColors = (/ 2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,24,26,28,0,37,39,41,43,45,47,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65 /)
  ;map@cnExplicitLabelBarLabelsOn = True 
  ;map@lbLabelStrings             =  (/"-1.00","-0.95","-0.90","-0.85","-0.80","-0.75","-0.70","-0.65","-0.60","-0.55","-0.50","-0.45","-0.40","-0.35","-0.30","-0.25","-0.20","-0.15","-0.10","-0.05","0.05","0.10","0.15","0.20","0.25","0.30","0.35","0.40","0.45","0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90","0.95","1.00"/)

;map@cnLevelSelectionMode = "ExplicitLevels"
  ;map@cnLevels     = (/ -10,-9.5,-9,-8.5,-8,-7.5,-7,-6.5,-6,-5.5,-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10 /)
  ;map@cnExplicitLabelBarLabelsOn = False 
  ;map@cnFillOpacityF = 1.0
  ;map@cnFillColors = (/ 2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,24,26,28,0,37,39,41,43,45,47,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65 /)
  ;map@cnExplicitLabelBarLabelsOn = True 
  ;map@lbLabelStrings             = (/"-10","-9.5","-9.0","-8.5","-8.0","-7.5","-7.0","-6.5","-6.0","-5.5","-5.0","-4.5","-4.0","-3.5","-3.0","-2.5","-2.0","-1.5","-1.0","-0.5","0.5","1.0","1.5","2.0","2.5","3.0","3.5","4.0","4.5","5.0","5.5","6.0","6.5","7.0","7.5","8.0","8.5","9.0","9.5","10"/)

;map@cnLevelSelectionMode = "ExplicitLevels"
  ;map@cnLevels     = (/  -2,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0 /)
  ;map@cnFillColors = (/ 2,3,3,4,4,5,6,7,8,9,0,14,15,16,17,18,19,19,20,20,21 /)
  ;map@cnExplicitLabelBarLabelsOn = True 
  ;map@lbLabelStrings             = (/ "-2.0","-1.8","-1.6","-1.4","-1.2","-1.0","-0.8","-0.6","-0.4","-0.2","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0" /)

;map@cnLevelSelectionMode = "ExplicitLevels"
  ;map@cnLevels     = (/ -0.003,-0.0027,-0.024,-0.0021,-0.0018,-0.0015,-0.0012,-0.0009,-0.0006,-0.0003,0.00045,0.0006,0.0009,0.0012,0.0015,0.0018,0.0021,0.0024,0.0027,0.003 /)
  ;map@cnFillColors = (/ 2,3,3,4,4,5,6,7,8,9,0,14,15,16,17,18,19,19,20,20,21 /)
  ;map@cnExplicitLabelBarLabelsOn = True 
  ;map@lbLabelStrings             = (/"-30","-27","-24","-21","-18","-15","-12","-9","-6","-3","3","6","9","12","15","18","21","24","27","30"/)

;map@cnLevelSelectionMode = "ExplicitLevels"
  ;map@cnLevels     =  (/-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1 /)
  ;map@cnFillColors = (/ 40,46,52,58,64,70,76,82,88,94,100,106,0 /)

;----------------------------------------------------------------------;

;--------------------  Paramétrisation du relief  ---------------------;
plot = gsn_csm_contour(wks,RES(:,:),map)

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
;pres@gsFillColor      = "lightsalmon4"
pres@gsFillColor      = "Grey"
pres@gsLineThicknessF = 2.0

; dessine le relief
dum0 = gsn_add_polygon(wks,plot,xpoly,ypoly,pres)
;----------------------------------------------------------------------;

;----------------------------  Iso-Theta  -----------------------------;
filename = dossier+"T"+sprinti("%0.11i",step)
TH = asciiread(filename,(/NLEV,NXPT/),"double")
iso   = new(NLEV-1,graphic)
xline = new(NXPT,"double")
yline = new(NXPT,"double")

tres                  = True
tres@gsLineThicknessF = 1.0

; Dessine les iso-théta
do i = 0, NXPT-1
  xline(i) = x(NLEV-1,i) 
end do

th0_old = 10.d10
do j = 10, NLEV-2
  th0 = (TH(j,0)+TH(j-1,0))/2.d0
  if (th0 .lt. th0_old - 4.d0) then
    yline(0) = (z(j,0) + z(j-1,0))/2.d0
    do i = 1, NXPT-1
      ;alpha = (th0 - TH(j-1,i))/(TH(j,i) - TH(j-1,i))
      ;yline(i) = alpha*z(j,i) + (1-alpha)*z(j-1,i)
    end do
    th0_old = th0
    ;iso(j-1) = gsn_add_polyline(wks,plot,xline,yline,tres)
  end if
end do
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
