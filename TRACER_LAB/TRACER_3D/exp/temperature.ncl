;############   SCRIPT NCL POUR TRACER DES CHAMPS   #############;
;#                                                              #;
;# auteur : C.Colavolpe                                         #;
;# sujet  : Permet de visualiser les résultats                  #; 
;#          ${experience} via la commande suivante :            #;
;#      ncl 'xp="XP"''plan="P"' 'step=n' temperature.ncl        #;
;#                                                              #;
;################################################################;

;------------  Definition des varialbes  -------------;
dossier = xp+"/res/"+plan+"/"
nom     = "_"+plan

x_borne = 150
x_unite = 30
y_borne = 100
y_unite = 20
z_borne = 20
z_unite = 1
;-----------------------------------------------------;

;---  Téléchargement de données pour les graphes  ----;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;-----------------------------------------------------;

;----------  Efface les messages d'erreur  -----------;
err = NhlGetErrorObjectId()
setvalues err
  "errLevel" : "Fatal"
end setvalues
;-----------------------------------------------------;

begin

;---------------  Lecture des champs  ----------------;
NLEV  = numAsciiRow(dossier+"Z"+sprinti("%0.11i",0))
NXPT  = numAsciiCol(dossier+"Z"+sprinti("%0.11i",0))
z = asciiread(dossier+"Z"+sprinti("%0.11i",0),(/NLEV,NXPT/),"double")
x = asciiread(dossier+"X"+sprinti("%0.11i",0),(/NLEV,NXPT/),"double")
filename = dossier+"T"+sprinti("%0.11i",step)
DATA = asciiread(filename,(/NLEV,NXPT/),"double")
;print("")
;print("DIMMENSIONS DU DOMAINE :")
;print("NXPT = "+NXPT+" ; NLEV = "+NLEV)
;print("")
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
res@lbLabelFontHeightF         = 0.00925
res@pmLabelBarWidthF           = 0.7
res@pmLabelBarHeightF          = 0.03
;-----------------------------------------------------;

;---------------  Choix des couleurs  ----------------;
res@cnFillOn = True
wks  = gsn_open_wks("epsi","T"+step+nom)
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
  ;res@cnLevels     = (/-0.003,-0.0027,-0.0024,-0.0021,-0.0018,-0.0015,-0.0012,-0.0009,-0.0006,-0.0004,0.00045,0.0006,0.0009,0.0012,0.0015,0.0018,0.0021,0.0024,0.0027,0.003 /)
  ;res@cnExplicitLabelBarLabelsOn = True 
  ;res@lbLabelStrings             = (/"-3.0","-2.7","-2.4","-2.1","-1.8","-1.5","-1.2","-0.9","-0.6","-0.3","0.3","0.6","0.9","1.2","1.5","1.8","2.1","2.4","2.7","3.0"/)
   ;res@cnFillColors = (/ 2,8,16,24,32,40,48,56,64,72,0,168,176,184,192,200,208,216,224,232,240 /)
;-----------------------------------------------------;

;--------------  Création des graphes  ---------------;
;plot = gsn_csm_contour(wks,DATA(::-1,:),res)
plot = gsn_csm_contour(wks,DATA(:,:),res)
draw(plot)
frame(wks)
delete(wks)
;-----------------------------------------------------;

end
