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
;experience = "QLHY"
;experience = "QLHY_EE"
experience = "QLHY_QE"
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
map@tiXAxisString       = "x(km)"
map@tiXAxisFontHeightF  = 0.0125
map@tmXBMinorPerMajor    = 1
map@tmXBMode             = "Manual"
map@tmXBLabelFontAspectF = 1.5
map@tmXBLabelFontHeightF = .0085
map@trXMinF              = -160.
map@tmXBTickStartF       = -160.
map@tmXBTickSpacingF     = 20.
map@trXMaxF              = 160.
map@tmXBTickEndF         = 160.
map@tmXBMajorOutwardLengthF = 0.0
map@tmXBMinorOutwardLengthF = 0
map@tmXBMajorLengthF = 0.0035
map@tmXBMajorThicknessF = 0.5
map@tmXBMinorLengthF = 0

map@sfYArray             = z
map@vpHeightF            = 0.6
;map@tiYAxisSide         = "Left"
map@tiYAxisString       = "z(km)"
map@tiYAxisFontHeightF  = 0.0125
map@tmYLMinorPerMajor    = 1
map@tmYLMode             = "Manual"
map@tmYLLabelFontAspectF = 1.5
map@tmYLLabelFontHeightF = .0085
map@trYMinF              = 0.
map@tmYLTickStartF       = 0.
map@tmYLTickSpacingF     = 1
map@trYMaxF              = 14
map@tmYLTickEndF         = 14.
map@tmYLMajorOutwardLengthF = 0.0
map@tmYLMinorOutwardLengthF = 0.0
map@tmYLMajorLengthF = 0.0035
map@tmYLMajorThicknessF = 0.5
map@tmYLMinorLengthF = 0

map@tmBorderThicknessF = 0.8

; Légende
map@lbLabelBarOn               = True
map@lbOrientation              = "vertical" 
map@cnLineLabelBackgroundColor = "transparent"
map@lbBoxLinesOn               = False
map@lbLabelAutoStride          = False
map@cnLabelBarEndStyle         = "ExcludeOuterBoxes" 
map@lbLabelFontHeightF         =.00825
;map@lbTitleOn                 = True
;map@lbTitleString             = "x10~S~-3"
;map@lbTitleFontHeightF        = 0.025
;map@pmLabelBarOrthogonalPosF  = .001 
;map@pmLabelBarParallelPosF    = 0.51
map@pmLabelBarWidthF           = 0.04
map@pmLabelBarHeightF          = 0.625

; Marque les ligne de séparation sans afficher la valeur
map@cnLinesOn           = True
map@cnLineLabelPerimOn  = False
map@cnLabelMasking      = False
map@cnLevelFlags        = True
map@cnLineThicknessF    = 1
map@cnLineLabelsOn      = False
map@gsnContourNegLineDashPattern = 2
;----------------------------------------------------------------------;

;------------------------  Choix des couleurs  ------------------------;
map@cnFillOn = True
wks  = gsn_open_wks("epsi",var+sprinti("%0.11i",step))


; Général de [-max ; +max]
;gsn_define_colormap(wks,"temp_19lev")
;map@cnLevelSelectionMode = "ManualLevels"
;map@cnMaxLevelValF  = doubletointeger(10*max((/abs(max(RES(:,:))),abs(min(RES(:,:)))/)))/10.
;map@cnMinLevelValF  = -map@cnMaxLevelValF
;map@cnLevelSpacingF = (map@cnMaxLevelValF-map@cnMinLevelValF)/15
;map@cnFillColors = (/ 2,3,4,5,6,7,8,9,0,14,15,16,17,18,19,20,21 /)

map@cnLevelSelectionMode = "ExplicitLevels"
gsn_define_colormap(wks,"MPL_Greys")
  map@cnLevels     = (/ -0.04,-0.035,-0.03,-0.025,-0.02,-0.015,-0.01,-0.005,0.,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04 /)
  map@cnExplicitLabelBarLabelsOn = False 
  map@cnFillOpacityF = 1.0
  map@cnFillColors = (/ 120,110,100,80,60,50,40,20,0,0,20,40,50,60,80,100,110 /)
  map@cnExplicitLabelBarLabelsOn = True 
  map@lbLabelStrings             = (/"-0.04", "-0.035","-0.03","-0.025","-0.02","-0.015","-0.01","-0.005","0.00","0.005","0.01","0.015","0.02","0.025","0.03","0.035", "0.04"/)
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
pres@gsFillColor      = "Black"
;pres@cnFillPattern     = 17
pres@gsLineThicknessF = 2.0

; dessine le relief
dum0 = gsn_add_polygon(wks,plot,xpoly,ypoly,pres)
;----------------------------------------------------------------------;

;----------------------------  Iso-Theta  -----------------------------;
;filename = dossier+"L"+sprinti("%0.11i",step)
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

;----------------------  Visionnage des graphes  ----------------------;
system("epstopdf "+var+sprinti("%0.11i",step)+".epsi")
system("evince "+var+sprinti("%0.11i",step)+".pdf")
system("rm -f *.epsi *.pdf")
system("ls -t "+dossier+var+"*")
;----------------------------------------------------------------------;

end
