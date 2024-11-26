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
experience = "COLD"
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
;map@tiXAxisString      = "Distance relative (en km)"
;map@tiXAxisFontHeightF = 0.018
map@vpWidthF           = 0.6
map@tmXBMinorPerMajor = 1
map@tmXBMode          = "Manual"
map@tmXBTickEndF      = 300
map@trXMaxF           = 300
map@tmXBTickEndF      = 25
map@trXMaxF           = 25
map@tmXBTickSpacingF  = 100
map@tmXBTickStartF    = -300
map@trXMinF           = -300


map@sfYArray = z
;map@tiYAxisSide        = "Left"
;map@tiYAxisString      = "Altitude (en km)"
;map@tiYAxisFontHeightF = 0.018
;map@vpHeightF          = 0.3
map@vpHeightF          = 0.3
map@tmYLMinorPerMajor = 1
map@tmYLMode          = "Manual"
;map@tmYLTickEndF      = 20
map@trYMaxF           = 21
map@tmYLTickEndF      = 21
map@tmYLTickSpacingF  = 1
map@tmYLTickStartF    = 1

; Légende
map@lbLabelBarOn               = True
map@lbOrientation              = "vertical" 
map@cnLineLabelBackgroundColor = "transparent"
map@lbBoxLinesOn               = True
map@lbLabelAutoStride          = False
map@cnLabelBarEndStyle         = "ExcludeOuterBoxes"
;map@lbTitleOn                  = True
;map@lbTitleString              = "x10~S~-3"
;map@lbTitleFontHeightF         = 0.025
;map@pmLabelBarOrthogonalPosF   = .001 
;map@pmLabelBarParallelPosF     = 0.51
map@pmLabelBarWidthF   = .04
map@pmLabelBarHeightF   = .345


; Marque les ligne de séparation sans afficher la valeur
map@cnLinesOn           = True
map@cnLineLabelPerimOn  = False
map@cnLabelMasking      = False
map@cnLevelFlags        = True
map@cnLineThicknessF    = 1.75
map@cnLineLabelsOn      = False
map@gsnContourNegLineDashPattern = 1
;----------------------------------------------------------------------;

;------------------------  Choix des couleurs  ------------------------;
map@cnFillOn = True
wks  = gsn_open_wks("epsi",var+sprinti("%0.11i",step))
gsn_define_colormap(wks,"temp_19lev")
;gsn_define_colormap(wks,"MPL_ocean")

; Général de [-max ; +max]
;map@cnLevelSelectionMode = "ManualLevels"
;  map@cnMaxLevelValF  = doubletointeger(10*max((/abs(max(RES(:,:))),abs(min(RES(:,:)))/)))/10.
;  map@cnMinLevelValF  = -map@cnMaxLevelValF
;  map@cnLevelSpacingF = (map@cnMaxLevelValF-map@cnMinLevelValF)/15
;  map@cnFillColors = (/ 2,3,4,5,6,7,8,9,0,14,15,16,17,18,19,20,21 /)

;map@cnLevelSelectionMode = "ExplicitLevels"
;  map@cnLevels     = (/  -2,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0 /)
;  map@cnFillColors = (/ 2,3,3,4,4,5,6,7,8,9,0,14,15,16,17,18,19,19,20,20,21 /)
;  map@cnExplicitLabelBarLabelsOn = True 
;  map@lbLabelStrings             = (/ "-2.0","-1.8","-1.6","-1.4","-1.2","-1.0","-0.8","-0.6","-0.4","-0.2","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0" /)

; HOT BUBLLE
map@cnLevelSelectionMode = "ExplicitLevels"
  map@cnLevels     = (/ -0.003,-0.00275,-0.0025,-0.00225,-0.0020,-0.00175,-0.0015,-0.00125,-0.001,-0.00075,-0.0005,0.0005,0.00075,0.001,0.00125,0.0015,0.00175,0.002,0.00225,0.0025,0.00275,0.003 /)
  map@cnFillColors = (/ 2,3,3,4,4,5,6,7,7,8,9,0,14,15,16,17,17,18,19,19,20,20,21 /)
  map@cnExplicitLabelBarLabelsOn = True 
  map@lbLabelStrings             = (/ "-3","-2.75","-2.5","-2.25","-2","-1.75","-1.5","-1.25","-1","-0.75 ","-0.5 ","0.5","0.75","1","1.25","1.5","1.75","2","2.25","2.5","2.75","3" /)

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
;----------------------------------------------------------------------;

;-----------------------  Création des graphes  -----------------------;
; dessine le relief
dum0 = gsn_add_polygon(wks,plot,xpoly,ypoly,pres)

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
