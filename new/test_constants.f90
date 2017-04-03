program p
use constants
implicit none
print*, "Eh_J       :", Eh_J       
print*, "Na         :", Na         
print*, "kcal_J     :", kcal_J     
print*, "Eh_kcalmol :", Eh_kcalmol 
print*, "a0_A       :", a0_A       
print*, "A_a0       :", A_a0       
print*, "c0         :", c0         
print*, "ea0        :", ea0        
print*, "Deb_ea0    :", Deb_ea0    
print*, "ea0_Deb    :", ea0_Deb    
print*, "Deb_eA     :", Deb_eA     
print*, "eA_Deb     :", eA_Deb     
print*, "h_Js       :", h_Js       
print*, "hbar_Js    :", hbar_Js    
print*, "Eh_cm1     :", Eh_cm1     
print*, "cm1_kcalmol:", cm1_kcalmol
print*, "cm1_eV     :", cm1_eV     
print*, "kB         :", kB         
print*, "e          :", e          
print*, "E_cc       :", E_cc       
print*, "CHARGECON  :", CHARGECON  
print*, "H_mass     :", H_mass     
print*, "O_mass     :", O_mass     

print*, ea0_Deb*A_a0
print*, eA_Deb
end program

