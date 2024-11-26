;#################  SCRIPT NCL POUR TRACER LES RMS ENTRE DEUX CHAMPS  ##################;
;#                                                                                     #;
;# auteurs : L.Auger, F.Voitus, C.Colavolpe, B.Menetrier                               #;
;# sujet   : Permet de visualiser les RMS des dossiers ${X1} et ${X2}                  #;
;#           via la commande :                                                         #;
;#                            ncl 'nom="exp"' 'var="X"' RMS.ncl                        #;
;#           où ${X3} une des variables {U,W,T,P}                                      #;
;#                                                                                     #;
;#######################################################################################;

;---------------------  Definition des varialbes  ---------------------;
nb_champs        = 100
NFREQOUTPUT_ARK  = 50
NFREQOUTPUT_UJ   = 100
NFREQOUTPUT_TRAP = 300
NFREQOUTPUT_RK4  = 300
RDT              = 0.1
;----------------------------------------------------------------------;

;------------  Téléchargement de données pour les graphes  ------------;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;----------------------------------------------------------------------;
begin
;----------------------  Definition des outils  -----------------------;
NLEV = numAsciiRow(nom+"/"+nom+"_RK4/res/Z"+sprinti("%0.11i",0))
NXPT = numAsciiCol(nom+"/"+nom+"_RK4/res/Z"+sprinti("%0.11i",0))
;----------------------------------------------------------------------;

;----------------  Valeurs des paramètres des graphes  ----------------;
res = True

res@gsnFrame = False
res@gsnDraw  = False

res@trYMinF = 0
res@trYMaxF = 0.1
res@trXMinF = 0
res@trXMaxF = nb_champs*0.1*NFREQOUTPUT_RK4
;----------------------------------------------------------------------;

;--------------------------  Calcul du RMS  ---------------------------;
do step = 2, nb_champs
  wks = gsn_open_wks("pdf","RMS_"+nom+"_"+sprinti("%0.3i",step)+".pdf")
  RMS = new((/3,step/),"double")
  X   = new(step,"double")

X(0)     = 0.d0
RMS(0,0) = 0.d0
RMS(1,0) = 0.d0
RMS(2,0) = 0.d0

do step_2 = 1, step-1

  X(step_2) = RDT*step_2*NFREQOUTPUT_RK4

  filename = nom+"/"+nom+"_RK4/res/"+var+sprinti("%0.11i",step_2*NFREQOUTPUT_RK4)
  ref  = asciiread(filename,(/NLEV,NXPT/),"double")

  filename = nom+"/"+nom+"_UJ/res/"+var+sprinti("%0.11i",step_2*NFREQOUTPUT_UJ)
  y = asciiread(filename,(/NLEV,NXPT/),"double")
  tmp = 0.d0
  do i = 0, NXPT-1
     do j = 0, NLEV-1
     	tmp = tmp + (y(j,i) - ref(j,i))^2
     end do
  end do
  RMS(0,step_2) = sqrt(tmp/NLEV/NXPT)

  filename = nom+"/"+nom+"_TRAP/res/"+var+sprinti("%0.11i",step_2*NFREQOUTPUT_TRAP)
  y = asciiread(filename,(/NLEV,NXPT/),"double")
  tmp = 0.d0
  do i = 0, NXPT-1
     do j = 0, NLEV-1
     	tmp = tmp + (y(j,i) - ref(j,i))^2
     end do
  end do
  RMS(1,step_2) = sqrt(tmp/NLEV/NXPT)

  filename = nom+"/"+nom+"_ARK/res/"+var+sprinti("%0.11i",step_2*NFREQOUTPUT_ARK)
  y = asciiread(filename,(/NLEV,NXPT/),"double")
  tmp = 0.d0
  do i = 0, NXPT-1
     do j = 0, NLEV-1
     	tmp = tmp + (y(j,i) - ref(j,i))^2
     end do
  end do
  RMS(2,step_2) = sqrt(tmp/NLEV/NXPT)

end do
;----------------------------------------------------------------------;
res@xyLineColors      = (/"blue","black","red"/)
res@xyLineThicknesses = (/3.0,3.0,3.0/)
res@xyDashPatterns   = (/0,0,0/)

xy = gsn_csm_xy(wks,X,RMS,res)
draw(xy)
frame(wks)
delete(wks)

delete(X)
delete(RMS)

  system("pdfcrop RMS_"+nom+"_"+ sprinti("%0.3i",step) + ".pdf RMS_"+nom+"_"+ sprinti("%0.3i",step)+".pdf")
end do
system("convert -delay 25 RMS_"+nom+"_* movie.gif")

system("eog movie.gif")
system("rm -f movie* *.pdf")
;----------------------------------------------------------------------;

end
