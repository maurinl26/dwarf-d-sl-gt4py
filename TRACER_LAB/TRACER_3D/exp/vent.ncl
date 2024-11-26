;#######   SCRIPT NCL POUR TRACER DES CHAMPS   ##########;
;#                                                      #;
;# auteur : C.Colavolpe                                 #;
;# sujet  : Permet de visualiser les résultats          #; 
;#          ${experience} via la commande suivante :    #;
;#      ncl 'xp="XP"' 'plan="P"' 'step=n' vent.ncl      #;
;#                                                      #;
;########################################################;

;------------  Definition des varialbes  -------------;
dossier    = xp+"/res/"+plan+"/"
nom        = "_"+plan

x_borne = 60
x_unite = 10
y_borne = 60
y_unite = 10
z_borne = 10
z_unite = 1
;-----------------------------------------------------;

;---  Téléchargement de données pour les graphes  ----;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;-----------------------------------------------------;

begin

;---------------  Lecture des champs  ----------------;
NLEV  = numAsciiRow(dossier+"Z"+sprinti("%0.11i",0))
NXPT  = numAsciiCol(dossier+"Z"+sprinti("%0.11i",0))
z = asciiread(dossier+"Z"+sprinti("%0.11i",0),(/NLEV,NXPT/),"double")
x = asciiread(dossier+"X"+sprinti("%0.11i",0),(/NLEV,NXPT/),"double")
filename = dossier+"U"+sprinti("%0.11i",step)
DATA = asciiread(filename,(/NLEV,NXPT/),"double")
print("")
print("DIMMENSIONS DU DOMAINE :")
print("NXPT = "+NXPT+" ; NLEV = "+NLEV)
print("")
;-----------------------------------------------------;

;-------------  Définition du résultat  --------------;
res = True

res@cnInfoLabelOn = False
res@gsnFrame = False
res@gsnDraw  = False
;-----------------------------------------------------;

;---------------  Définition des axes  ---------------;
res@cnFillMode = "AreaFill"
res@trGridType = "TriangularMesh"

; Abscisses
res@sfXArray           = x
res@vpWidthF           = 0.6
if (plan .eq. "Y") then
  res@tiXAxisString       = "~F14~y~F13~(km)"
else
  res@tiXAxisString       = "~F14~x~F13~(km)"
end if
res@tiXAxisFontHeightF = 0.015

res@tmXBMinorPerMajor       = 1
res@tmXBMode                = "Manual"
res@tmXBLabelFontAspectF    = 1.5
res@tmXBLabelFontHeightF    = 0.015
res@tmXBMajorOutwardLengthF = 0.005
res@tmXBMinorOutwardLengthF = 0
res@tmXBMajorLengthF        = 0.005
res@tmXBMinorLengthF        = 0
res@tmXBPrecision           = 2

res@trXMinF            = -x_borne
res@tmXBTickStartF     = -x_borne
res@tmXBTickSpacingF   =  x_unite
res@trXMaxF            =  x_borne
res@tmXBTickEndF       =  x_borne

; Ordonnées
res@sfYArray            = z
res@tiYAxisSide         = "Left"
if (plan .eq. "Z") then
  res@tiYAxisString       = "~F14~y~F13~(km)"
else
  res@tiYAxisString       = "~F14~z~F13~(km)"
end if
res@tiYAxisFontHeightF  = 0.015

res@tmYLMinorPerMajor       = 1
res@tmYLMode                = "Manual"
res@tmYLLabelFontAspectF    = 1.5
res@tmYLLabelFontHeightF    = 0.015
res@tmYLMajorOutwardLengthF = 0.005
res@tmYLMinorOutwardLengthF = 0
res@tmYLMajorLengthF        = 0.005
res@tmYLMinorLengthF        = 0
res@tmYLAutoPrecision       = False

if (plan .eq. "Z") then
  res@vpHeightF        = 0.6*y_borne/x_borne
  res@tmYLPrecision    = 2
  res@trYMinF          = -y_borne
  res@tmYLTickStartF   = -y_borne
  res@tmYLTickSpacingF =  y_unite
  res@trYMaxF          =  y_borne
  res@tmYLTickEndF     =  y_borne
else
  res@vpHeightF        = 0.2 + 0*0.6*z_borne/x_borne
  res@tmYLPrecision    = 2
  res@trYMinF          = 0
  res@tmYLTickStartF   = 0
  res@tmYLTickSpacingF = z_unite
  res@trYMaxF          = z_borne
  res@tmYLTickEndF     = z_borne
end if
;-----------------------------------------------------;

;---------------------  Légende  ---------------------;
res@lbLabelBarOn               = True
res@lbOrientation              = "horizontal" 
res@cnLineLabelBackgroundColor = "transparent"
res@lbBoxLinesOn               = False
res@lbLabelAutoStride          = False
res@cnLabelBarEndStyle         = "ExcludeOuterBoxes" 
res@lbLabelFontHeightF         =.00925
res@pmLabelBarWidthF           = 0.7
res@pmLabelBarHeightF          = 0.03
;-----------------------------------------------------;

;---------------  Choix des couleurs  ----------------;
res@cnFillOn = True
wks  = gsn_open_wks("epsi","U"+step+nom)
gsn_define_colormap(wks,"BlueWhiteOrangeRed")

res@cnLinesOn          = True
res@cnLineLabelPerimOn = False
res@cnLabelMasking     = True
res@cnLevelFlags       = True
res@cnLineThicknessF   = 2
res@cnLineLabelsOn     = False
res@cnFillOpacityF     = 1.0
res@gsnContourNegLineDashPattern = 1

; Général de [-max ; +max]
res@cnLevelSelectionMode = "ManualLevels"
  res@cnMaxLevelValF  = doubletointeger(10*max((/abs(max(DATA(:,:))),abs(min(DATA(:,:)))/)))/10.
  res@cnMinLevelValF  = -res@cnMaxLevelValF
  res@cnLevelSpacingF = (res@cnMaxLevelValF-res@cnMinLevelValF)/15
  res@cnFillColors = (/ 2,12,22,32,42,52,62,72,0,165,174,183,192,201,210,219,228 /)

;res@cnLevelSelectionMode = "ExplicitLevels"
  ;res@cnLevels     = (/-1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1 /)
  ;res@cnLevels     = (/  -2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0 /)
  ;res@cnLevels     = (/ -3,-2.85,-2.7,-2.55,-2.4,-2.25,-2.1,-1.95,-1.8,-1.65,-1.5,-1.35,-1.2,-1.05,-0.9,-0.75,-0.6,-0.45,-0.3,-0.15,0.15,0.3,0.45,0.6,0.75,0.9,1.05,1.2,1.35,1.5,1.65,1.8,1.95,2.1,2.25,2.4,2.55,2.7,2.85,3 /)
  ;res@cnExplicitLabelBarLabelsOn = False 
  ;res@cnFillColors = (/ 2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,24,26,28,0,37,39,41,43,45,47,49,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65 /)
  ;res@cnExplicitLabelBarLabelsOn = True 
  ;res@lbLabelStrings             = (/"-1.00","-0.95","-0.90","-0.85","-0.80","-0.75","-0.70","-0.65","-0.60","-0.55","-0.50","-0.45","-0.40","-0.35","-0.30","-0.25","-0.20","-0.15","-0.10","-0.05","0.05","0.10","0.15","0.20","0.25","0.30","0.35","0.40","0.45","0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90","0.95","1.00"/)
  ;res@lbLabelStrings             = (/ "-2.0","-1.9","-1.8","-1.7","-1.6","-1.5","-1.4","-1.3","-1.2","-1.1","-1.0","-0.9","-0.8","-0.7","-0.6","-0.5","-0.4","-0.3","-0.2","-0.1","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0","1.1","1.2","1.3","1.4","1.5","1.6","1.7","1.8","1.9","2.0" /)
;-----------------------------------------------------;

;------------  Paramétrisation du relief  ------------;
;plot = gsn_csm_contour(wks,DATA(::-1,:),res)
plot = gsn_csm_contour(wks,DATA(:,:),res)

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
pres@gsFillColor      = "Grey"
pres@gsLineThicknessF = 2.0

; dessine le relief
;dum0 = gsn_add_polygon(wks,plot,xpoly,ypoly,pres)
;-----------------------------------------------------;

;--------------------  Iso-Theta  --------------------;
;filename = dossier+"T"+sprinti("%0.11i",step)
;TH = asciiread(filename,(/NLEV,NXPT/),"double")
;iso   = new(NLEV-1,graphic)
;xline = new(NXPT,"double")
;yline = new(NXPT,"double")

;tres                  = True
;tres@gsLineThicknessF = 1.0

; Dessine les iso-théta
;do i = 0, NXPT-1
;  xline(i) = x(NLEV-1,i) 
;end do

;th0_old = 10.d10
;do j = 10, NLEV-2
;  th0 = (TH(j,0)+TH(j-1,0))/2.d0
;  if (th0 .lt. th0_old - 4.d0) then
;    yline(0) = (z(j,0) + z(j-1,0))/2.d0
;    do i = 1, NXPT-1
;      alpha = (th0 - TH(j-1,i))/(TH(j,i) - TH(j-1,i))
;      yline(i) = alpha*z(j,i) + (1-alpha)*z(j-1,i)
;    end do
;    th0_old = th0
;    iso(j-1) = gsn_add_polyline(wks,plot,xline,yline,tres)
;  end if
;end do
;----------------------------------------------------------------------;

;-----------------------  Création des graphes  -----------------------;
draw(plot)
frame(wks)
delete(wks)
;----------------------------------------------------------------------;

end
