;###############  SCRIPT NCL POUR VISION DYNAMIQUE DES CHAMPS CALCULÉS  ################;
;#                                                                                     #;
;# auteurs : L.Auger, F.Voitus, C.Colavolpe, B.Menetrier                               #;
;# sujet   : Permet de visualiser les résultats du dossier ${experience} de            #;
;#           l'état initial à ${NTIMESTEP} via la commande suivante :                  #;
;#                               ncl 'var="X"' animate.ncl                             #;
;#           où ${X} une des variables {U,W,T,PI,P}                                    #;
;#                                                                                     #;
;#######################################################################################;

;---------------------  Definition des varialbes  ---------------------;
experience = "NLNH/NLNH_RK4"
NTIMESTEP   = 300*100
NFREQOUTPUT = 300
RDT         = 0.1
dossier = "./"+experience+"/res/"
;dossier = "../res/"
;----------------------------------------------------------------------;

;------------  Téléchargement de données pour les graphes  ------------;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;----------------------------------------------------------------------;

;-----------------------  Dessin des champs  --------------------------;
begin

REF   = NTIMESTEP
NLEV  = numAsciiRow(dossier+"X"+sprinti("%0.11i",0))
NXPT  = numAsciiCol(dossier+"X"+sprinti("%0.11i",0))

x = asciiread(dossier+"X"+sprinti("%0.11i",0),(/NLEV,NXPT/),"double")

mov   = NTIMESTEP/NFREQOUTPUT
tmp   = new((/NLEV,NXPT/),"double")
field = new((/mov+1,NLEV,NXPT/),"double")

filename = dossier+var+sprinti("%0.11i",REF)
tmp(:,:) = asciiread(filename,(/NLEV,NXPT/),"double")
max_f    = max(tmp(:,:))
min_f    = min(tmp(:,:))

;----------------  Valeurs des paramètres des graphes  ----------------;
map = True

; Enlève des infos inutiles
map@cnInfoLabelOn = False

map@gsnDraw  = False
map@gsnFrame = False

; Axes abscisse et ordonnée
map@cnFillMode = "AreaFill"
map@trGridType = "TriangularMesh"
map@sfXArray   = x
  ;map@tiXAxisString      = "Distance relative (en km)"
  ;map@tiXAxisFontHeightF = 0.018
  map@vpWidthF           = 0.6

z  = asciiread(dossier+"Z"+sprinti("%0.11i",REF),(/NLEV,NXPT/),"double")
map@sfYArray = z
  ;map@tiYAxisSide        = "Left"
  ;map@tiYAxisString      = "Altitude (en km)"
  ;map@tiYAxisFontHeightF = 0.018
  map@vpHeightF          = 0.6

  map@tmXBMinorPerMajor = 1
  map@tmXBMode          = "Manual"
  map@tmXBTickEndF      = 500.
  map@tmXBTickSpacingF  = "50"
  map@tmXBTickStartF    = "-100."
  map@tmYLMinorPerMajor = 1
  map@tmYLMode          = "Manual"
  map@tmYLTickEndF      = 18
  map@tmYLTickSpacingF  = "1"
  map@tmYLTickStartF    = "1"

; Légende
map@lbLabelBarOn               = False
map@lbOrientation              = "vertical" 
map@cnLineLabelBackgroundColor = "transparent"
map@lbBoxLinesOn               = False
  map@lbLabelAutoStride     = False 

; Marque les ligne de séparation sans afficher la valeur
map@cnLinesOn           = True
map@cnLineLabelPerimOn  = False
map@cnLabelMasking      = False
map@cnLevelFlags        = True
map@cnLineThicknessF    = 2.3
map@cnLineLabelsOn      = False
map@cnRasterSmoothingOn = False
map@gsnContourNegLineDashPattern = 0
;----------------------------------------------------------------------;

;------------------------  Choix des couleurs  ------------------------;
map@cnFillOn = True

; Général de [-max ; +max]
;map@cnLevelSelectionMode = "ManualLevels"
;map@cnMaxLevelValF  = doubletointeger(10*max((/abs(max(RES(:,:))),abs(min(RES(:,:)))/)))/10.
;map@cnMinLevelValF  = -map@cnMaxLevelValF
;map@cnLevelSpacingF = (map@cnMaxLevelValF-map@cnMinLevelValF)/15
;map@cnFillColors = (/ 2,3,4,5,6,7,8,9,0,14,15,16,17,18,19,20,21 /)

map@cnLevelSelectionMode = "ExplicitLevels"
  map@cnLevels     = (/  -2,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0 /)
  map@cnFillColors = (/ 2,3,3,4,4,5,6,7,8,9,0,14,15,16,17,18,19,19,20,20,21 /)

map@cnExplicitLabelBarLabelsOn = True 
map@lbLabelStrings             = (/ "-2.0","-1.8","-1.6","-1.4","-1.2","-1.0","-0.8","-0.6","-0.4","-0.2","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0" /)
;----------------------------------------------------------------------;

;------------------------  Création des films  ------------------------;
do step=0,mov
  wks = gsn_open_wks("pdf","movie_"+ var + "_" + sprinti("%0.3i",step)+".pdf")

  gsn_define_colormap(wks,"ncl_default")

  filename        = dossier+var+sprinti("%0.11i",step*NFREQOUTPUT)
  tmp(:,:)        = asciiread(filename,(/NLEV,NXPT/),"double")
  field(step,:,:) = tmp(:,:)
  z  = asciiread(dossier+"Z"+sprinti("%0.11i",step*NFREQOUTPUT),(/NLEV,NXPT/),"double")
  map@sfYArray = z
  heures   = floattointeger(step*NFREQOUTPUT*RDT/3600)
  minutes  = floattointeger((step*NFREQOUTPUT*RDT-3600*heures)/60)
  secondes = floattointeger(step*NFREQOUTPUT*RDT-3600*heures-60*minutes)
  if (heures .lt. 10) then
    heure = "0"+integertoshort(heures)
  else
    heure = integertoshort(heures)
  end if
  if (minutes .lt. 10) then
    minute = "0"+integertoshort(minutes)
  else
    minute = integertoshort(minutes)
  end if
  if (secondes .lt. 10) then
    seconde = "0"+integertoshort(secondes)
  else
    seconde = integertoshort(secondes)
  end if
  ;map@tiMainString = "Variable "+var+" apres "+heure+"h "+minute+"min "+seconde+"s"
  plot = gsn_csm_contour(wks,field(step,:,:),map)

  draw(plot)
  frame(wks)
  delete(wks)

  system("pdfcrop movie_"+ var + "_" + sprinti("%0.3i",step) + ".pdf movie_"+ var + "_" + sprinti("%0.3i",step)+".pdf")
end do
;----------------------------------------------------------------------;

;----------------------  Visionnage des graphes  ----------------------;
;system("rm -f *.pdf")
system("convert -delay 25 movie_"+var+"* movie_"+var+".gif")

system("eog movie_"+var+".gif")
system("rm -f movie_*")
;----------------------------------------------------------------------;

end
