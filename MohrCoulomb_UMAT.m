function [STRESS,hsv,DDSDDE]=MohrCoulomb_UMAT(pseudo_time,PROPS,STRESS0,DSTRAN0,hsv)
% [sig_new,hsv0,DDSDDE]=DP(pseudo_time,Props,sig,d_strain,hsv);
warning off;
format long;
format compact;
%% parameters
eta    = 500;
E      =  PROPS(1);   %!  E    - Young's modulus
nu     =  PROPS(2);   %!  nu   - Poisson's ratio
phi   =  PROPS(3);  %!  phi  - angle of friction (rad)
% phi_d   = c_reduction(phi_d,10,eta,hsv);
psi   =  PROPS(4);  %!  psi  - angle of dilation (rad)
coh     =  PROPS(5);  %!  c    - cohesion
% coh   = c_reduction(coh,0.1,eta,hsv);
%% change the order of stress and strain 
nsigma=4;
STRESS=STRESS0(1:4);
DSTRAN=zeros(nsigma,1);

%
DSTRAN(1:2)=DSTRAN0(1:2);
DSTRAN(4)=DSTRAN0(3);
DSTRAN(3)=0;
% %% initialization
% if pseudo_time==0
%     hsv=zeros(50,1);
% end
%% Elastic stiffness

[D,Dinv]=DlinElas(E,nu,nsigma);

DSTRES = D*DSTRAN; %! stress increment
SigB   = STRESS + DSTRES;  %! elastic predictor stress
% !-----------------------------------------------------------------------
% !     Stress update and creation of consistent constitutive matrix
% !-----------------------------------------------------------------------
PlastPar(1) = (1.0+sin(phi))/(1.0-sin(phi));  %! Friction parameter
PlastPar(2) = 2.0*coh*sqrt(PlastPar(1));        %! Uniaxial comp. str
PlastPar(3) = (1.0+sin(psi))/(1.0-sin(psi));  %! Dilatation parameter
%% return mapping
[SigC,Depc,region]=MohrCoulombStressReturn(SigB,nsigma,PlastPar,D,Dinv);
%% calculate the plastic strain
% delta_sig=SigB-SigC;
% stain_els=D\delta_sig;
% deps_plas=DSTRAN-stain_els;
DelStrainPlas = (Dinv*(SigB-SigC));%! Plastic strain increment
DPE_eq = sqrt(dot(DelStrainPlas,DelStrainPlas)*0.6666666666666667); %! Equivalent plastic strain increment
% !-----------------------------------------------------------------------
% !     Postconditioner for ABAQUS
% !-----------------------------------------------------------------------
% DDSDDE     = Depc;          %! consistent constitutive matrix

STRESS = [SigC; 0; 0];
% !-----------------------------------------------------------------------
% !	Type of failure and Lode angle
% !-----------------------------------------------------------------------
[depsp_v, depsp_d]=strain_invariant(DelStrainPlas,nsigma); % DelStrainPlas shear strain component NOT multiply by 2
% hsv(1)  = region;	          % ! type of stress return
% hsv(17) = hsv(17)+depsp_v;    % ! Equivalent plastic strain increment
hsv = hsv+depsp_d;     % ! Total equivalent plastic strain
%% return the DDSDDE consistent constitutive matrix
DDSDDE(1,1)=Depc(1,1);
DDSDDE(2,2)=Depc(2,2);
% DDSDDE(4,4)=Depc(3,3);
DDSDDE(3,3)=Depc(4,4);
%
DDSDDE(1,2)=Depc(1,2);
DDSDDE(2,1)=DDSDDE(1,2);
% DDSDDE(1,4)=Depc(1,3);
% DDSDDE(4,1)=DDSDDE(1,4);
%
% DDSDDE(2,4)=Depc(2,3);
% DDSDDE(4,2)=DDSDDE(2,4);
%%

end
%%
function [Sigma_up,Depc,region]=MohrCoulombStressReturn(Sigma,nsigma,PlasPar,D,Dinv)
% !------------------------------------------------------------------------
% ! Wrapping subroutine for return mapping with a linearly elastic perfectly
% ! plastic Mohr-Coulomb material. Performed in principal coordinates.
% !
% ! The Mohr-Coulomb criterion is defined by
% !		f = k*sigP(1) - sigP(3) - comp = 0
% !
% ! and the corresponding plastic potential
% !       g = m*sigP(1) - sigP(3)
% !
% !
% ! INPUT
% !  - Name -		-type and description --
% !	Sigma		real array(nsigma). Stress vector. Size is dependent on the type of stress.
% !				  The ordering of the components must be as follows:
% !				    Plane situations (plane stress, plane strain and axisymmetry):
% !					Sigma = [sig_x sig_y sig_z tau_xy] (i.e. sig_z is out-of-plane)
% !					General 3D:  Sigma = [sig_x sig_y sig_z tau_xy tau_xz tau_yz]
% !	nsigma		integer,scalar. Number of stress components, i.e. the size of Sigma
% !	PlasPar		real, array(3). Vector of plastic material parameters [k comp m]
% !							k: Describes internal friction, k = (1+sin(phi))/(1-sin(phi))
% !							comp: Uniaxial compressive strength. comp = 2*c*sqrt(k)
% !							m: Describes dilation. m = (1+sin(psi))/(1-sin(psi))
% !	D			real, array(nsigma,nsigma) Elastic constitutive matrix.
% !	Dinv		real, array(nsigma,nsigma) Inverted elastic constitutive matrix (compliance matrix).
% !
% ! OUTPUT
% !	Sigma_up	real, array(nsigma) Stress vector containg the updated stresses.
% !	Depc		real, array(nsigma,nsigma) Consistent constitutive elasto-plastic matrix
% !				  If the material is not yielding Depc = D
% !	region		integer, scalar. Number of the region of the returned stress
% !					region = 0: the elastic region, i.e. no yielding
% !					region = 1: Return to the Mohr-Coulomb surface
% !					region = 2: Return to the triaxial compressive line, sigp_1 = sigp_2
% !					region = 3: Return to the triaxial tensile line, sigp_2 = sigp_3
% !					region = 4: Return to the apex, sigp_1 = sigp_2 = sigp_3.
% !
% !------------------------------------------------------------------------
% ! Johan Clausen
% ! Department of Civil Engineering
% ! Aalborg University
% ! May 2008
% !------------------------------------------------------------------------
% 	implicit none
% 	integer(4), intent(in)  :: nsigma
% 	real(8), intent(in)     :: Sigma(nsigma), PlasPar(3)
% 	real(8), intent(in)     :: D(nsigma,nsigma), Dinv(nsigma,nsigma)
%
% 	real(8), intent(out)	:: Sigma_up(nsigma), Depc(nsigma,nsigma)
% 	integer(4), intent(out) :: region
% 	!-----------------------
%
% 	integer(4) hoop,ouplP, s1, s2, nshear
% 	real(8) psi(3,3), f, SigP(3), SigP_up(3)
% 	real(8), dimension(nsigma,nsigma) :: DepP, DepcP, A, Atrans
% 	real(8), dimension(nsigma,nsigma) :: MatProd, T
% 	real(8) Fnorm(3), Gnorm(3), Lfdir(3), Lgdir(3), SiPla(3)
% 	real(8) Tshear(nsigma-3,nsigma-3),Tshearinv(nsigma-3,nsigma-3)

% 	!--- Settings regarding calculation of the constitutive matrix, DepP -------------
DepFormLine = 1;
DepFormApex = 1;
beta = 20.0;
alpha = 100.0;
% 	integer(4), parameter :: DepFormLine = 1 ! Form of the consistent constitutive matrix on a line
% 					! DepFormLine = 0: Double singular constitutive matrix
% 					! DepFormLine = 1: Modified "almost double singular" matrix on a line ("Tiln?_b=axr" in the MatLab code)
% 					! DepFormLine = 2: A single singular matrix in the Koiter direction, i.e. the direction of the plastic strain.
% 	real(8), parameter :: beta = 20.0_8 ! The value used in the modification of DepP when DepFormLine = 1
%
% 	integer(4), parameter :: DepFormApex = 1 ! Form of the consistent constitutive matrix on the apex
% 					! DepFormApex = 0: DepP is the zero matrix
% 					! DepFormApex = 1: A modified single singular DepP in the Koiter direction ("DepKo/fak" in the MatLab code)
% 					! DepFormApex = 2: A modified double singular DepP in the Koiter and the hydrostatic direction ("DepKoiteP" in the MatLab code)
% 	real(8), parameter :: alpha = 100.0_8 ! ! Factor in the modified apex-Depc-formulation. Used when DepFormApex = 1 or DepFormApex = 2
% 	!--------------------------------------------------------------------------------
%
% 	!---- Explicit interfaces ----------------------
% 	interface
% 		pure subroutine PrinStressAna(Sigma,nsigma,SigP)
% 			integer(4), intent(in) :: nsigma
% 			real(8), intent(in) :: Sigma(nsigma)
% 			real(8), intent(out) :: SigP(3)
% 		end subroutine PrinStressAna
%
% 		pure subroutine PrinRetMoCo(SigP,f,PlasPar,Dfull,nsigma,
%      1                SigP_up,region)
% 			integer(4), intent(in) :: nsigma
% 			real(8), intent(in)    :: SigP(3), PlasPar(3), f
% 			real(8), intent(in)    :: Dfull(nsigma,nsigma)
% 			real(8), intent(out)  :: SigP_up(3)
% 			integer(4), intent(out) :: region
% 		end subroutine PrinRetMoCo
%
% 		pure subroutine TshearPrinPerfect(SigP,SigP_up,nshear,s1,s2,
%      1												Tshear,Tshearinv)
% 			real(8), intent(in) :: SigP(3), SigP_up(3)
% 			integer(4), intent(in) :: nshear, s1, s2
% 			real(8), dimension(nshear,nshear), intent(out) :: Tshear, Tshearinv
% 		end subroutine TshearPrinPerfect
%
% 		pure subroutine FormDepPerfect(D,Norm,Edir,nsigma,Dep)
% 			integer(4), intent(in) :: nsigma
% 			real(8), intent(in) :: D(nsigma,nsigma), Norm(nsigma)
% 			real(8), intent(in) :: Edir(nsigma)
% 			real(8), intent(out) :: Dep(nsigma,nsigma)
% 		end subroutine FormDepPerfect
%
% 		pure subroutine FormModDepLine(mod_type,beta,SiPla,Fline,Gline,Dninv,
%      1											Dc,Dcinv,nsigma,Depc)
% 			integer(4), intent(in) :: mod_type, nsigma
% 			real(8), intent(in) :: Dninv(3,3),Dc(nsigma,nsigma),Gline(3)
% 			real(8), intent(in) :: beta,Sipla(3),Fline(3),Dcinv(nsigma,nsigma)
% 			real(8), intent(out) :: Depc(nsigma,nsigma)
% 		end subroutine FormModDepLine
%
% 		pure subroutine FormModDepApex(mod_type,alpha,SiPla,Dninv,
%      &										Dc,Dcinv,nsigma,Depc)
% 			integer(4), intent(in) :: mod_type, nsigma
% 			real(8), intent(in) :: alpha,Sipla(3),Dninv(3,3)
% 			real(8), intent(in) :: Dc(nsigma,nsigma),Dcinv(nsigma,nsigma)
% 			real(8), intent(out) :: Depc(nsigma,nsigma)
% 		end subroutine FormModDepApex
%
% 		pure subroutine PrinDirect(Sigma,nsigma,SigP,psi)
% 			integer(4), intent(in) :: nsigma
% 			real(8), intent(in) :: Sigma(nsigma), SigP(3)
% 			real(8), intent(out)  :: psi(3,3)
% 		end subroutine PrinDirect
%
% 		pure subroutine TransMatrix(psi,nsigma,ouplP,A)
% 			integer(4), intent(in) :: nsigma, ouplP
% 			real(8), intent(in)    :: psi(3,3)
% 			real(8), intent(out)   :: A(nsigma,nsigma)
% 		end subroutine TransMatrix
%
% 	end interface
% 	!------------------------------------------
% !------------------------------------------------------------------------
[SigP]=PrinStressAna(Sigma,nsigma); %! Principal stresses SigP

% !	--- Value of Yield function f ---
% !	PlasPar = [k comp m]
f = PlasPar(1)*SigP(1) - SigP(3) - PlasPar(2); %! f = k*SigP_1 - SigP_3 - comp
if (f > 0.0)
    % !		--- Position of out-of plane or hoop stress in principal stress vector SigP ----
    if (nsigma == 4)
        if (Sigma(3) == SigP(1))
            ouplP = 1 ;
            s1 = 2 ;
            s2 = 3;
        elseif (Sigma(3) == SigP(2))
            ouplP = 2 ;
            s1 = 1 ;
            s2 = 3;
        elseif (Sigma(3) == SigP(3))
            ouplP = 3 ;
            s1 = 1 ;
            s2 = 2;
        end
    end
    % !		---------------------------------------------------------------------------------
    % !		--- Stress return in principal coordinates --
    [SigP_up,region]=PrinRetMoCo(SigP,f,PlasPar,D,nsigma);
    SiPla = SigP - SigP_up; %! Plastic corrector stress in principal stress space
    % !		---------------------------------------------
    % !		--- Modification matrix, T ----------------------
    nshear = nsigma - 3; %! Number of shear components
    [Tshear,Tshearinv]=TshearPrinPerfect(SigP,SigP_up,nshear,s1,s2); %! Forms the shear part modification matrix T
    T = zeros(nsigma,nsigma);
    T(1:3,1:3) = eye(3); %reshape(1.0, 0.0, 0.0,0.0, 1.0, 0.0, 0.0, 0.0, 1.0); %! Unit Matrix
    T(4:nsigma,4:nsigma) = Tshear;
    % !		-------------------------------------------------
    % !		--- Infinitesimal constitutive matrix, DepP -----------------
    DepP = 0.0;
    if (region == 1)   % ! Return to the yield surface
        DepP(4:nsigma,4:nsigma) = D(4:nsigma,4:nsigma);	  %! Shear components of consistent constitutive matrix in principal stress space
        
        % 			! PlasPar = [k comp m]
        Fnorm(1:3) = ([PlasPar(1) 0.0 -1.0]); %! Yield plane normal
        Gnorm = ([PlasPar(3), 0.0, -1.0]);    %! Potential plane normal
        
        [DepP(1:3,1:3)]=FormDepPerfect(D(1:3,1:3),Fnorm,Gnorm,3); %! Normal components of consistent constitutive matrix in principal stress space
        
    elseif (region == 2)  %! return to a line 1, triaxial compression, sigp1 = sigp2
        Lfdir = ([1.0, 1.0, PlasPar(1)]); %! Edge line direction
        Lgdir = ([1.0, 1.0, PlasPar(3)]); %! Potential edge line direction
        
        [DepP]=FormModDepLine(DepFormLine,beta,SiPla,Lfdir,Lgdir,Dinv(1:3,1:3),D,Dinv,nsigma);
        
    elseif (region == 3 )  %! return to a line 1, triaxial extension, sigp2 = sigp3
        Lfdir = ([1.0, PlasPar(1), PlasPar(1)]); %! Edge line direction
        Lgdir = ([1.0, PlasPar(3), PlasPar(3)]); %! Potential edge line direction
        
        [DepP]=FormModDepLine(DepFormLine,beta,SiPla,Lfdir,Lgdir,Dinv(1:3,1:3),D,Dinv,nsigma);
        
    else %! apex return
        [DepP]=FormModDepApex(DepFormApex,alpha,SiPla,Dinv(1:3,1:3),D,Dinv,nsigma);
    end  %! region = ...
    % !		------------------------------------------------------------
    
    % !		--- Consistent constitutive matrix, DepcP -------
    DepcP = T*DepP;
    % !		-------------------------------------------------
    
    % !		--- Tranformation matrix A ------------------
    [psi]=PrinDirect(Sigma,nsigma,SigP);
    [A]=TransMatrix(psi,nsigma,ouplP);
    
    % !		--- Coordinate transformation ---------------
    Atrans = (A');
    Sigma_up = zeros(nsigma,1);
    for i=1:nsigma
        Sigma_up(i,1) =dot(SigP_up,Atrans(i,1:3));
    end
    MatProd = (Atrans*DepcP);
    Depc = (MatProd*A);
    % !		---------------------------------------------
    
else %! no yielding
    Sigma_up = Sigma;
    Depc = D;
    region = 0;
end  %! f > 0

% !--------------------------------------------------------------------------
end %subroutine MohrCoulombStressReturn !-----------------------------
% !--------------------------------------------------------------------------
%%
function [SigP_up,region]=PrinRetMoCo(SigP,f,PlasPar,Dfull,nsigma)
% !-----------------------------------------------------------------------
% ! Stress return for Mohr-Coulomb plasticity in principal stress space.
% ! Includes non-associated plasticity, but no hardening is allowed.
% ! The Mohr-Coulomb criterion is written as f = k*sigP_1 - sigP_3 - comp = 0
% ! The plastic potential is m = m*sigP_1 - sigP_3
% !
% ! INPUT
% !  - Name - -type,size -  -- Description --
% !     SigP  (real,ar(3))  Principal predictor stresses in descending order
% !                   SigP = [sigP_1 sigP_2 sigP_3]
% !     f	  (real,sc)   Value of the yield function given by
% !                   f = k*sigP_1 - sigP_3 - comp
% !     PlasPar (real,ar(3))  Vector of plastic material parameters  = [k comp m]
% !                   k: Describes internal friction, k = (1+sin(phi))/(1-sin(phi))
% !                   comp: Uniaxial compressive strength. comp = 2*c*sqrt(k)
% !                   m: Describes dilation. m = (1+sin(psi))/(1-sin(psi))
% !                   psi = Friction angle, c = cohesion, psi = dilation angle
% !
% !     Dfull (real,ar(nsigma,nsigma)) Elastic constitutive matrix
% !     nsigma  (integer,sc)  Number of stress components, Dfull
% !
% ! OUTPUT
% !     SigP_up (real,ar(3))  Updated principal stresses, i.e. they obey the
% !                   Mohr-Coulomb yield criterion
% !     region  (integer,sc)  Type of return, i.e. the region of the
% !                   predictor stress
% !
% !----------------------------------------------------------------------
% ! Johan Clausen
% ! Esbjerg Technical Institute
% ! Aalborg University
% ! November 2005
% !----------------------------------------------------------------------

% %       !-------------------
%       integer(4), intent(in) :: nsigma
%       real(8), intent(in)    :: SigP(3), PlasPar(3), f
%       real(8), intent(in)    :: Dfull(nsigma,nsigma)
%
%       real(8), intent(out)  :: SigP_up(3)
%       integer(4), intent(out) :: region
% %       !-------------------
%
%       integer(4) s1, s2
%       real(8) D(3,3)
%       real(8) k, m, comp, den, Rp(3), RpHat(2)
%       real(8) Rp2(3), Rp3(3), N2(3), N3(3), t1, t2
%       real(8) NI_II(3), NI_III(3), pI_II, pI_III
%       real(8) apex, SigPApex(3), den_t
% !------------------------------------------------------------------------

%       !--- Preliminary parameters ----------------
k  = PlasPar(1);   %! Friction parameter
m  = PlasPar(3);   %! Dilation parameter
comp = PlasPar(2); %! Uniaxial compressive strength
apex = comp/(k-1); %! Stress coordinate of the criterions apex
if (isinf(apex))>0
    apex=comp;
end
D = Dfull(1:3,1:3);%! Relates to normal stresses


den = k*(D(1,1)*m-D(1,3)) - D(3,1)*m + D(3,3); %! denominator a'*D*b
Rp(1) = (D(1,1)*m-D(1,3))/den; %! Rp is the scaled direction of
Rp(2) = (D(2,1)*m-D(2,3))/den; %! the update direction,
Rp(3) = (D(3,1)*m-D(3,3))/den; %! Rp = D*b/(a'*D*b) a, b gradient of
%                       ! yield surface and plastic potential, respectively.
% !-------------------------------------------

%       !--- Boundary planes -----------------------
SigPapex = SigP - apex; %! Vector from predictor stress to the apex

%       ! Boundary plane between regions I and II
NI_II(1) = Rp(2)*k - Rp(3);   %! NI_II = cross(Rp,R1)
NI_II(2) = Rp(3)   - Rp(1)*k; %! R1 = [1 1 k], direction
NI_II(3) = Rp(1)   - Rp(2);   %! vector of line 1

pI_II = NI_II(1)*SigPapex(1) + NI_II(2)*SigPapex(2)+ NI_II(3)*SigPapex(3);

%       ! Boundary plane between regions I and III
NI_III(1) = Rp(2)*k - Rp(3)*k; %! NI_III = cross(Rp,R2)
NI_III(2) = Rp(3)   - Rp(1)*k; %! R2 = [1 k k], direction
NI_III(3) = Rp(1)*k - Rp(2);   %! vector of line 2

pI_III = NI_III(1)*SigPapex(1) + NI_III(2)*SigPapex(2)+ NI_III(3)*SigPapex(3);
%       !-------------------------------------------

%       !--- t-paramters for region determination --
%         ! secondary surface in region II a = [0 k -1], b = [0 m -1]
den = k*(D(2,2)*m-D(2,3)) - D(3,2)*m + D(3,3); %! denominator a'*D*b
Rp2(1) = (D(1,2)*m-D(1,3))/den; %! Rp is the scaled direction of
Rp2(2) = (D(2,2)*m-D(2,3))/den; %! the update direction,
Rp2(3) = (D(3,2)*m-D(3,3))/den; %! Rp = D*b/(a'*D*b) a, b gradient of
%                        ! yield surface and plastic potential, respectively.

N2(1) = Rp(2)*Rp2(3) - Rp(3)*Rp2(2); %! N2 = cross(Rp,Rp2)
N2(2) = Rp(3)*Rp2(1) - Rp(1)*Rp2(3);
N2(3) = Rp(1)*Rp2(2) - Rp(2)*Rp2(1);

den_t = N2(1) + N2(2) + k*N2(3); %! N2'*R1, R1 = [1 1 k] Direction of line 1
%       ! t-parameter of line 1
t1 = (N2(1)*SigPapex(1) + N2(2)*SigPapex(2)+ N2(3)*SigPapex(3))/den_t;


%       ! secondary surface in region III a = [k -1 0], b = [m -1 0]
den = k*(D(1,1)*m-D(1,2)) - D(2,1)*m + D(2,2); %! denominator a'*D*b
Rp3(1) = (D(1,1)*m-D(1,2))/den; %! Rp is the scaled direction of
Rp3(2) = (D(2,1)*m-D(2,2))/den; %! the update direction,
Rp3(3) = (D(3,1)*m-D(3,2))/den; %! Rp = D*b/(a'*D*b) a, b gradient of
%                        ! yield surface and plastic potential, respectively.

N3(1) = Rp(2)*Rp3(3) - Rp(3)*Rp3(2); %! N3 = cross(Rp,Rp3)
N3(2) = Rp(3)*Rp3(1) - Rp(1)*Rp3(3);
N3(3) = Rp(1)*Rp3(2) - Rp(2)*Rp3(1);

den_t = N3(1) + k*N3(2) + k*N3(3); %! N3'*R2, R2 = [1 k k] Direction of line 2
%       ! t-parameter of line 1
t2 = (N3(1)*(SigP(1)-apex) + N3(2)*(SigP(2)-apex)+ N3(3)*(SigP(3)-apex))/den_t;
%       !-------------------------------------------

%       !--- Region determination and update -------
if (t1 > 0.0 && t2 > 0.0)
    region = 4;
    SigP_up(1:3) = apex;
elseif (pI_II < 0.0)
    region = 2;
    SigP_up(1) = t1   + apex; %! SigP_up = t1*R1 + apex
    SigP_up(2) = t1   + apex; %! R1 = [1 1 k], direction
    SigP_up(3) = t1*k + apex; %! vector of line 1
elseif (pI_III <= 0.0)
    region = 1;
    SigP_up(1) = SigP(1) - f*Rp(1); %! SigP_up = SigP - SiPla
    SigP_up(2) = SigP(2) - f*Rp(2); %! SiPla is the plastic corrector
    SigP_up(3) = SigP(3) - f*Rp(3); %! given by f*Rp
else
    region = 3;
    SigP_up(1) = t2   + apex; %! SigP_up = t2*R2 + apex
    SigP_up(2) = t2*k + apex; %! R2 = [1 k k], direction
    SigP_up(3) = t2*k + apex; %! vector of line 2
end  %! Regions and update
%       !----------------------------------------------

% !--------------------------------------------------------------------------------
end %subroutine PrinRetMoCo !-----------------------------------------------
% !--------------------------------------------------------------------------------
%%
function Depc=FormModDepLine(mod_type,beta,SiPla,Fline,Gline,Dninv,Dc,Dcinv,nsigma)
% !--------------------------------------------------------------------------
% ! Calculates the modified elasto-plastic constitutive matrix at on a line
% ! for a perfectly plastic material, It is modified in the sence that is not
% ! the double singular matrix that theory prescribes. Only one type of
% ! modification is possible. The calculation are carried out in principal
% ! stress space.
% !
% ! INPUT
% ! - Name -		-type and description --
% !	mod_type	(integer,scalar) Detemines the modified type
% !				  mod_type= 0: A Doubly singular matrix, as given in the
% !							   References [1] and [2].
% !				  mod_type= 1: A single singular matrix in the plastic
% !							   strain direction, but almost singular in
% !							   any direction perpendicular to the potential
% !							   line. See Reference [3] for a detailed explanation.
% !	beta		(real,scalar) The factor by which to divide the stiffness in
% !				  the directions perpendicular to the line.
% !	SiPla		(real,array(3)) The plastic corrector principal stresses.
% !	Fline		(real,array(3)) Direction of the yield surface line in principal
% !				  stress space.
% !	Gline		(real,array(3)) Direction of the plastic potential line in principal
% !				  stress space.
% !	Dninv		(real,arrey(3,3)) Normal components of the elastic
% !				  compliance matrix.
% !	Dc			(real,array(nsigma,nsigma)) Modified elastic stiffness.
% !				  If the yield criterion is linear Dc = D.
% !	Dcinv		(real,array(nsigma,nsigma)) Inverse of the modified
% !				  elastic stiffness. If the yield criterion is linear,
% !				  Dcinv = Dinv.
% !	nsigma		(integer, scalar) Number of stress components
% !
% ! OUTPUT
% !	Depc	(real,array(nsigma,nsigma) Modified Elasto-plastic constitutive
% !			  matrix on a line. If the criterion is linear Depc is the
% !			  infinitesimal version, otherwise it is the consistent. In
% !			  both cases it is expressed in principal stres space.
% !
% ! REFERENCES:
% ! [1] Johan Clausen, Lars Damkilde and Lars Andersen: "Efficient return
% !	  algorithms for associated plasticity with multiple yield planes",
% !	  International Journal for Numerical Methods in Engineering, 2006,
% !	  volume 66, pages 1036-1059.
% ! [2] Johan Clausen, Lars Damkilde and Lars Andersen: "An efficient return
% !	  algorithm for non-associated plasticity with linear yield criteria
% !	  in principal stress space". Computers & Structures, 2007, volume 85,
% !	  pages 1795-1807.
% ! [3] Johan Clausen: Efficient Non-Linear Finite Element Implementation of
% !	  Elasto-Plasticity for Geotechnical Problems. Ph.D. Thesis. 2006.
% !	  Esbjerg Technical Institute, Aalborg University. Can be downloaded
% !	  from http://vbn.aau.dk/fbspretrieve/14058639/JCthesis.pdf
% !
% !------------------------------------------------------------------------
% ! Johan Clausen
% ! Department of Civil Engineering
% ! Aalborg University
% ! April 2008
% !------------------------------------------------------------------------
% 	implicit none
% 	!----- Input/Output ------
% 	integer(4), intent(in) :: mod_type, nsigma
% 	real(8), intent(in) :: Dninv(3,3),Dc(nsigma,nsigma),Gline(3)
% 	real(8), intent(in) :: beta,Sipla(3),Fline(3),Dcinv(nsigma,nsigma)
%
% 	real(8), intent(out) :: Depc(nsigma,nsigma)
% 	!-------------------------
%
% 	real(8) KoitDir(3),KoitPerDir(3), Dper(3,3)
%
% 	!---- Explicit interfaces ----------------------
% 	interface
% 		pure subroutine FormDepPerfect(D,Norm,Edir,nsigma,Dep)
% 			integer(4), intent(in) :: nsigma
% 			real(8), intent(in) :: D(nsigma,nsigma), Norm(nsigma)
% 			real(8), intent(in) :: Edir(nsigma)
% 			real(8), intent(out) :: Dep(nsigma,nsigma)
% 		end subroutine FormDepPerfect
% 		!------
% 		pure subroutine Cross(A,B,C)
% 			real(8), intent(in)	 ::	A(3), B(3)
% 			real(8), intent(out) :: C(3)
% 		end subroutine Cross
% 		!-----
% 		pure subroutine FormDepLinePerfect(D,Dinv,Ra,Rb,nsigma,DepLine)
% 			integer(4), intent(in) :: nsigma
% 			real(8), intent(in), dimension(nsigma,nsigma) :: D, Dinv
% 			real(8), intent(in), dimension(3) :: Ra, Rb
% 			real(8), intent(out) :: DepLine(nsigma,nsigma)
% 		end subroutine FormDepLinePerfect
% 	end interface
% 	!-----------------------------------------------
% !------------------------------------------------------------------------
Depc = zeros(nsigma,nsigma);

Depc=FormDepLinePerfect(Dc,Dcinv,Fline,Gline,nsigma);

if (mod_type == 1 && dot(SiPla,SiPla) > 0.0)
    KoitDir = (SiPla*Dninv); %! the plastic strain direction
    KoitPerDir=Cross(KoitDir,Gline); %! Direction perpendicular to the strain direction
    % 											 ! and the plastic potential line
    Dper=FormDepLinePerfect(Dc(1:3,1:3),Dcinv(1:3,1:3),KoitPerDir,KoitPerDir,3);
    Depc(1:3,1:3) = Depc(1:3,1:3) + Dper/beta;
elseif (mod_type == 2 && dot(SiPla,SiPla) > 0.0)
    KoitDir = (SiPla*Dninv); %! The plastic strain direction in principal stress space
    Depc(1:3,1:3)=FormDepPerfect(Dc(1:3,1:3),KoitDir,KoitDir,3);
    Depc(4:nsigma,4:nsigma) = Dc(4:nsigma,4:nsigma);
end

% !--------------------------------------------------------------------------------
end %subroutine FormModDepLine !--------------------------------------------
% !--------------------------------------------------------------------------------
%%
function Depc=FormModDepApex(mod_type,alpha,SiPla,Dninv,Dc,Dcinv,nsigma)
% !--------------------------------------------------------------------------
% ! Calculates the modified elasto-plastic constitutive matrix at the apex
% ! for a perfectly plastic material, It is modified in the sence that is not
% ! the zero matrix. Two different modifications are possible, and they are
% ! determined by the variable mod_type. The calculations are carried out in
% ! principal stress space.
% !
% ! INPUT
% ! - Name -		-type and description --
% !	mod_type	(integer,scalar) Detemines the modified type
% !				  mod_type= 0: Depc is the zero matrix.
% !				  mod_type= 1: A single singular matrix in the plastic
% !							   strain direction.
% !				  mod_type= 2: A double singular matrix in the plastic
% !							   strain direction and in the direction of
% !							   the hydrostatic stress line. This should
% !							   be used for apex points not on the hydrostatic
% !							   line (ishydro = 0)
% !	alpha		(real,scalar) The facto by which to divide the stiffness
% !	SiPla		(real,array(3)) The plastic corrector principal stresses.
% !	Dninv		(real,arrey(3,3)) Normal components of the elastic
% !				  compliance matrix.
% !	Dc			(real,array(nsigma,nsigma)) Modified elastic stiffness.
% !				  If the yield criterion is linear Dc = D.
% !	Dcinv		(real,array(nsigma,nsigma)) Inverse of the modified
% !				  elastic stiffness. If the yield criterion is linear,
% !				  Dcinv = Dinv.
% !	nsigma		(integer, scalar) Number of stress components
% !
% ! OUTPUT
% !	Depc	(real,array(nsigma,nsigma) Modified Elasto-plastic constitutive
% !			  matrix on an apex. If the criterion is linear Depc is the
% !			  infinitesimal version, otherwise it is the consistent. In
% !			  both cases it is expressed in principal stres space.
% !
% !------------------------------------------------------------------------
% ! Johan Clausen
% ! Departmen of Civil Engineering
% ! Aalborg University
% ! April 2008
% !------------------------------------------------------------------------
% implicit none
% integer(4), intent(in) :: mod_type, nsigma
% real(8), intent(in) :: alpha,Sipla(3),Dninv(3,3)
% real(8), intent(in) :: Dc(nsigma,nsigma),Dcinv(nsigma,nsigma)
%
% real(8), intent(out) :: Depc(nsigma,nsigma)
% !-------------------------
%
% real(8) KoitDir(3), Dnorm(3,3), Pdir(3), PerDir(3)
%
% !---- Explicit interfaces ----------------------
% interface
% pure subroutine FormDepPerfect(D,Norm,Edir,nsigma,Dep)
% integer(4), intent(in) :: nsigma
% real(8), intent(in) :: D(nsigma,nsigma), Norm(nsigma)
% real(8), intent(in) :: Edir(nsigma)
% real(8), intent(out) :: Dep(nsigma,nsigma)
% end subroutine FormDepPerfect
% !--------
% pure subroutine FormDepLinePerfect(D,Dinv,Ra,Rb,nsigma,DepLine)
% integer(4), intent(in) :: nsigma
% real(8), intent(in), dimension(nsigma,nsigma) :: D, Dinv
% real(8), intent(in), dimension(3) :: Ra, Rb
% real(8), intent(out) :: DepLine(nsigma,nsigma)
% end subroutine FormDepLinePerfect
% !------
% pure subroutine Cross(A,B,C)
% real(8), intent(in)	 ::	A(3), B(3)
% real(8), intent(out) :: C(3)
% end subroutine Cross
% end interface
% !------------------------------------------
% !-----------------------------------------------------------------------
Depc = zeros(nsigma,nsigma);

if (dot(SiPla,SiPla) > 0.0)
    KoitDir = SiPla*Dninv;
    
    %     !if (mod_type == 0) then the Zero matrix is used
    if (mod_type == 1)  %! Single singular matrix (DepKo/fak)
        Dnorm=FormDepPerfect(Dc(1:3,1:3),KoitDir,KoitDir,3);
        Depc(1:3,1:3) = Dnorm/alpha; %! Normal part
        Depc(4:nsigma,4:nsigma) = Dc(4:nsigma,4:nsigma)/alpha; %! Shear part
    elseif (mod_type == 2)  %! Double singular matrix (KoiteP)
        Pdir = ([1.0, 1.0, 1.0]); %! Hydrostatic direction
        PerDir=Cross(KoitDir,Pdir);
        Depc=FormDepLinePerfect(Dc,Dcinv,Perdir,PerDir,nsigma);
        Depc = Depc/alpha;
    end
end  %! dot_product(SiPla,SiPla) > 0.0_8

%         !--%------------------------------------------------------------------------------
end % subroutine FormModDepApex !--------------------------------------------
%     !------%--------------------------------------------------------------------------
%%
function [Tshear,Tshearinv]=TshearPrinPerfect(SigP,SigP_up,nshear,s1,s2)
% !--------------------------------------------------------------------------
% ! Forms the shear part of the modification matrix T in principal stress space
% ! for at perfectly plastic material. T is used when calculation the modified
% ! elastic stiffness matrix Dc which is, in turn, used to form the consistent
% ! constitutive matrix. The shear part of T is unaffected whether the yield
% ! criterion is linear or not.
% !
% !
% ! INPUT
% !  - Name -	-type,size -	-- Description --
% !	SigP	 real,ar(3)		Principal predictor stresses in descending order
% !							  SigP = [sigP_1 sigP_2 sigP_3]
% !	SigP_up	 real,ar(3)		Updated principal stresses in descending order
% !	nshear	 int,sc			Number of shear stress components. 1 in plane problems
% !							  and 3 in full 3D
% !	s1		 int,sc			Position of largest principal in-plane stress in SigP.
% !							  Only used in plane situations
% !	s2		 int,sc			Position of smallest principal in-plane stress in SigP.
% !							  Only used in plane situations
% !
% ! OUTPUT
% !	Tshear		real,ar(nshear,nshear) Shear part of the modification matrix T
% !	Tshearinv	real,ar(nshear,nshear) Inverse of Tshear
% !
% !------------------------------------------------------------------------
% ! Johan Clausen
% ! Department of Civil Engineering
% ! Aalborg University
% ! April 2008
% !------------------------------------------------------------------------
% 	implicit none
% 	real(8), intent(in) :: SigP(3), SigP_up(3)
% 	integer(4), intent(in) :: nshear, s1, s2
%
% 	real(8), dimension(nshear,nshear), intent(out) :: Tshear, Tshearinv
% !------------------------------------------------------------------------

Tshear = 0.0;
Tshearinv = 0.0;

if (nshear == 1)
    if ((SigP_up(s1) - SigP_up(s2)) > 0.0)
        Tshear = (SigP_up(s1)-SigP_up(s2)) / (SigP(s1)-SigP(s2));
        Tshearinv = 1.0/Tshear;
    end
elseif (nshear == 3)
    if ((SigP_up(1) - SigP_up(2)) > 0.0  && (SigP(1) - SigP(2)) > 0.0)
        Tshear(1,1) = (SigP_up(1)-SigP_up(2)) / (SigP(1)-SigP(2));
        Tshearinv(1,1) = 1.0/Tshear(1,1);
    end
    if ((SigP_up(1) - SigP_up(3)) > 0.0  && (SigP(1) - SigP(3)) > 0.0)
        Tshear(2,2) = (SigP_up(1)-SigP_up(3)) / (SigP(1)-SigP(3));
        Tshearinv(2,2) = 1.0/Tshear(2,2);
    end
    if ((SigP_up(2) - SigP_up(3)) > 0.0  &&   (SigP(2) - SigP(3)) > 0.0)
        Tshear(3,3) = (SigP_up(2)-SigP_up(3)) / (SigP(2)-SigP(3));
        Tshearinv(3,3) = 1.0/Tshear(3,3);
    end
end

% !--------------------------------------------------------------------------------
end %subroutine TshearPrinPerfect !-----------------------------------------
% !--------------------------------------------------------------------------------
%%
function SigP=PrinStressAna(Sigma,nsigma)
% !--------------------------------------------------------------------------
% ! Calculates the principal stress of a stress vector Sigma. Analytical
% ! expressions are used in the calculations.
% !
% ! INPUT
% !  - Name - -type,size -  -- description --
% !     Sigma (real,ar(nsigma)Stress vector. Size is dependent on the type of stress
% !                   nsigma = 4 in plane situations and nsigma = 6 in 3D
% !                   The ordering of the components must be as follows:
% !                   Plane situations (plane stress, plane strain and axisymmetry):
% !                   Sigma = [sig_x sig_y sig_z tau_xy] (i.e. sig_z is out-of-plane)
% !                   General 3D:  Sigma = [sig_x sig_y sig_z tau_xy tau_xz tau_yz]
% !     nsigma  (integer,sc)  Number of stress components, i.e. the size of Sigma
% !
% ! OUTPUT
% !     SigP  (real,ar(3))  Principal stresses in descending order
% !
% !---------------------------------------------------------------------------
% ! Johan Clausen
% ! Esbjerg Institute of Technology
% ! Aalborg University
% ! November 2005
% !--------------------------------------------------------------------------
% 	implicit none
%       !----- Input/Output -----------
%       integer(4), intent(in) :: nsigma
%       real(8), intent(in) :: Sigma(nsigma)
%
%       real(8), intent(out) :: SigP(3)
% 	!------------------------------
%
%       real(8) sig_av, sig_hj, I1, J2, J3, lode, sigm
%       real(8) sqJ2, sin3lode
%       real(8) S(6)
%
%       !---- Explicit interfaces ----------------------
% 	interface
% 		pure subroutine Invariants(Sigma,nsigma,S,I1,J2,J3,lode,sin3lode)
% 			integer(4), intent(in) :: nsigma
% 			real(8), intent(in) :: Sigma(nsigma)
% 			real(8), intent(out) :: S(nsigma), I1, J2
% 			real(8), intent(out) :: J3, lode, sin3lode
% 		end subroutine Invariants
% 	end interface
% 	!-----------------------------------------------
%
% !------------------------------------------------------------------------
% SigP = zeros(3,1);

% !----- Plane situations including axisymmetry -----------------------------

if (nsigma == 4)
    sig_av = 0.5*(Sigma(1) + Sigma(2)) ;
    sig_hj = sqrt((0.5*(Sigma(1) - Sigma(2)))^2 + Sigma(4)^2);
    SigP(1) = sig_av + sig_hj;
    SigP(2) = sig_av - sig_hj;
    SigP(3) = Sigma(3); %! The out-of-plane stress
    %       !-- Sorting: SigP(1) >= SigP(2) >= SigP(3) --
    if (Sigma(3) > SigP(1))
        SigP(3) = SigP(2)  ;
        SigP(2) = SigP(1)  ;
        SigP(1) = Sigma(3);
    elseif (Sigma(3) > SigP(2) )
        SigP(3) = SigP(2)  ;
        SigP(2) = Sigma(3) ;
    end
    % !----- General stress state ---------------------------------------------
elseif (nsigma == 6)
    % !   ----- Invariants ----------------------------------------------------
    [S,I1,J2,J3,lode,sin3lode]=Invariants(Sigma,nsigma);
    sigm = I1/3.0; %! Hydrostatic stress
    
    if (sqrt(J2) < 1.0e-12*abs(I1))
        SigP(1) = Sigma(1);
        SigP(2) = Sigma(2);
        SigP(3) = Sigma(3);
    else
        % !   ----- Principal stresses -----------------------------------------
        sqJ2 = 1.154700538379252*sqrt(J2);
        SigP(1) = sqJ2 * dsin(lode + 2.094395102393195) + sigm;
        SigP(2) = sqJ2 * dsin(lode)                     + sigm;
        SigP(3) = sqJ2 * dsin(lode - 2.094395102393195) + sigm;
    end
    
end  %! nsigma == ...
% !----------------------------------------------------------------------------------
end %subroutine PrinStressAna !-----------------------------------------------
% !----------------------------------------------------------------------------------
%%
function Dep=FormDepPerfect(D,Norm,Edir,nsigma)
% !--------------------------------------------------------------------------
% ! Calculates elasto-plastic constitutive matrix of a material with
% ! non-associated and perfect plasticity. If Norm = Edir the flow rule
% ! is associated
% !
% ! INPUT
% !     D   (real) elastic constitutive matrix. Size = nsigma x nsigma
% !     Norm  (real) Gradient (normal) of the yield surface. Size = nsigma
% !     Edir  (real) Gradient (normal) of the plastic potential. Size = nsigma
% !     nsigma  (integer) Number of stress components
% !
% ! OUTPUT
% !     Dep   (real) Elasto-plastic infinitesimal constitutive matrix
% !
% !--------------------------------------------------------------------------
% ! Johan Clausen
% ! Esbjerg Department of Engineering
% ! Aalborg University
% ! July 2005
% !--------------------------------------------------------------------------
% 	implicit none
%       integer(4), intent(in) :: nsigma
% 	real(8), intent(in) :: D(nsigma,nsigma), Norm(nsigma)
% 	real(8), intent(in) :: Edir(nsigma)
%       real(8), intent(out) :: Dep(nsigma,nsigma)
% % !--------------------------------------------------------------------------
%       real(8), dimension(nsigma) :: Num1, Num2
%       real(8) den, Num(nsigma,nsigma)
%       integer(4) i,j
% !--------------------------------------------------------------------------
Num1 = Edir*D;
Num2 = Norm*D;
for i = 1:nsigma
    for j = 1:nsigma
        Num(i,j) = Num1(i)*Num2(j);
    end
end
den  = dot(Norm,Num1);
Dep = D - Num/den;
% !--------------------------------------------------------------------------
end %subroutine FormDepPerfect !--------------------------------------
% !--------------------------------------------------------------------------
%%
function DepLine=FormDepLinePerfect(D,Dinv,Ra,Rb,nsigma)
% !--------------------------------------------------------------------------
% ! Calculates double singular elasto-plastic constitutive matrix of a
% ! material with non-associated and perfect plasticity. If Ra = Rb the flow
% ! rule is associated. The formula is only valid in principal stress space
% !
% ! INPUT
% !     D   (real) elastic constitutive matrix. Size = nsigma x nsigma
% !     Dinv  (real) Inverted elastic constitutive matrix. Size = nsigma x nsigma
% !     Ra    (real) Direction of the line on the yield surface in principal
% !           stress space. If the line is curved Ra is the tangent
% !     Rb    (real) Direction of the line on the plastic potential in principal
% !           stress space. If the line is curved Rb is the tangent
% !
% ! OUTPUT
% !     Depline (real) Double singular elasto-plastic infinitesimal
% !           constitutive matrix on a line. Size = nsigma x nsigma
% !
% !------------------------------------------------------------------------
% ! Johan Clausen
% ! Esbjerg Department of Engineering
% ! Aalborg University
% ! July 2005
% !------------------------------------------------------------------------

%       integer(4), intent(in) :: nsigma
%       real(8), intent(in), dimension(nsigma,nsigma) :: D, Dinv
%       real(8), intent(in), dimension(3) :: Ra, Rb
%
%       real(8), intent(out) :: DepLine(nsigma,nsigma)
% %       !-------------------------
%
%       real(8) den, Den1(3), Num(3,3), Dprin(3,3)
%       integer(4) i,j
for i = 1:3
    for j = 1:3
        Num(i,j) = Ra(i)*Rb(j);
    end
end

Den1 = Rb*Dinv(1:3,1:3);
den = dot(Ra,Den1);        %! denomimator, Ra'*Dinv*Rb

Dprin = Num/den;

DepLine = zeros(nsigma,nsigma);
DepLine(1:3,1:3) = Dprin;
if (nsigma > 3)
    DepLine(4,4) = D(4,4); %! In plane situations nsigma == 4
    if (nsigma == 6)       %! General three dimensional stress state
        DepLine(5,5) = D(5,5);
        DepLine(6,6) = D(6,6);
    end
end
% !-------------------------------------------------------------------------------
end %subroutine FormDepLinePerfect !---------------------------------------
% !-------------------------------------------------------------------------------
%%
function psi=PrinDirect(Sigma,nsigma,SigP)
% !--------------------------------------------------------------------------
% ! Calculates the principal directions of a principal stress state
% ! Analytical expressions are used in the calculations.
% !
% ! INPUT
% !  - Name -	-type,size -	-- description --
% !	Sigma	(real,ar(nsigma)Stress vector. Size is dependent on the type of stress
% !							  nsigma = 4 in plane situations and nsigma = 6 in 3D
% !							  The ordering of the components must be as follows:
% !							  Plane situations (plane stress, plane strain and axisymmetry):
% !							  Sigma = [sig_x sig_y sig_z tau_xy] (i.e. sig_z is out-of-plane)
% !							  General 3D:  Sigma = [sig_x sig_y sig_z tau_xy tau_xz tau_yz]
% !	nsigma	(integer,sc)	Number of stress components, i.e. the size of Sigma
% !	SigP	(real,ar(3))	Principal stresses in descending order
% !
% ! OUTPUT
% !	psi		(real,ar(3,3))	Principal angle if nsigma = 4 or
% !							  vector of normalized eigenvectors (principal directions)
% !							  if nsigma = 6. When nsigma = 4 the angle is stored
% !							  in the upper left corner, i.e. psi(1,1).
% !
% !---------------------------------------------------------------------------
% ! Johan Clausen
% ! Esbjerg Institute of Technology
% ! Aalborg University
% ! November 2005
% !--------------------------------------------------------------------------
% 	integer(4), intent(in) :: nsigma
% 	real(8), intent(in) :: Sigma(nsigma), SigP(3)
%
% 	real(8), intent(out)  :: psi(3,3)
%
% !----- Internal variables -------------------------------------------------
% 	real(8), parameter :: pi = 3.14159265358979323_8
% 	real(8) tol1, tol2, tol, nl(3), normal(3,3), leng
% 	real(8) SigMax, SigMin
% 	integer(4) i, hj1, hj2, nr
%
% %       !---- Explicit interfaces ----------------------
% 	interface
% 		pure subroutine Cross(A,B,C)
% 			real(8), intent(in)	 ::	A(3), B(3)
% 			real(8), intent(out) :: C(3)
% 		end subroutine Cross
% 	end interface
% 	!-----------------------------------------------
% !----- Plane situations including axisymmetry -----------------------------

if (nsigma == 4)
    
    if (Sigma(1) > Sigma(2)  &&  Sigma(4) >= 0.0)
        psi = 0.5*atan(2*Sigma(4)/(Sigma(1)-Sigma(2)));
    elseif (Sigma(1) < Sigma(2)  &&  Sigma(4) >= 0.0)
        psi = 0.5*(pi - atan(2*Sigma(4)/(Sigma(2)-Sigma(1))));
    elseif (Sigma(1) < Sigma(2)  &&  Sigma(4) < 0.0)
        psi = 0.5*(atan(-2*Sigma(4)/(Sigma(2)-Sigma(1))) + pi);
    elseif (Sigma(1) > Sigma(2)  &&  Sigma(4) < 0.0)
        psi = 0.5*(2*pi - atan(-2*Sigma(4)/(Sigma(1)-Sigma(2))));
    elseif (Sigma(1) == Sigma(2) && Sigma(4) > 0.0)
        psi = 0.25*pi;
    elseif (Sigma(1) == Sigma(2) && Sigma(4) < 0.0)
        psi = 0.75*pi;
    elseif (Sigma(1) == Sigma(2) && Sigma(4) == 0.0)
        psi = 0.0 ;
    end
    
    % !------------------------------------------------------------------------
    
    
    % !----- General 3D ----------------------------------------------------
elseif (nsigma == 6)
    tol1 = abs(1.0e-10*(SigP(1) - SigP(3))); %! tolerance value 1
    tol2 = abs(1.0e-12*(SigP(1) + SigP(2) + SigP(3))); %! tolerance value 2
    tol = max(([tol1,tol2,1.0e-13]));
    
    nl = 0.0;
    normal = 0.0;
    
    % !	-- Sigma is already expressed in principal coordinates ---
    if (abs(Sigma(4)) < tol  &&  abs(Sigma(5)) < tol  && abs(Sigma(6)) < tol)
        if (Sigma(1) >= Sigma(2) && Sigma(1) >= Sigma(3))
            normal(1,1) = 1.0;
            if (Sigma(2) >= Sigma(3))
                normal(2,2) = 1.0 ;
                normal(3,3) = 1.0;
            else
                normal(3,2) = 1.0 ;
                normal(2,3) = 1.0;
            end
        elseif (Sigma(2) >= Sigma(1) && Sigma(2) >= Sigma(3))
            normal(2,1) = 1.0;
            if (Sigma(1) >= Sigma(3))
                normal(1,2) = 1.0 ;
                normal(3,3) = 1.0;
            else
                normal(3,2) = 1.0 ;
                normal(1,3) = 1.0;
            end
        elseif (Sigma(3) >= Sigma(1) && Sigma(3) >= Sigma(2))
            normal(3,1) = 1.0;
            if (Sigma(1) >= Sigma(2))
                normal(1,2) = 1.0 ;
                normal(2,3) = 1.0;
            else
                normal(2,2) = 1.0 ;
                normal(1,3) = 1.0;
            end
        end
        % !	-- Three distinct eigenvalues -----
    elseif (abs(SigP(1)-SigP(2)) > tol  && abs(SigP(2)-SigP(3)) > tol  && abs(SigP(1)-SigP(3)) > tol)
        for i = 1:2
            nl(1) = (Sigma(2)-SigP(i))* (Sigma(3)-SigP(i)) - Sigma(6)^2;
            nl(2) = -Sigma(4) * (Sigma(3)-SigP(i))+ Sigma(6) * Sigma(5);
            nl(3) = Sigma(4) * Sigma(6)- (Sigma(2)-SigP(i)) * Sigma(5);
            leng = sqrt(nl(1)^2 + nl(2)^2 + nl(3)^2);
            if (leng < tol)
                normal(i,i) = 1 ;
            else
                normal(:,i) = nl/leng;
            end
        end
        % !			normal(:,3) = cross(Normal(:,1),Normal(:,2))
        normal(:,3)=cross(normal(:,1),normal(:,2));
        % !	-- Three equal eigenvalues ----
    elseif (abs(SigP(1)-SigP(2)) < tol  && abs(SigP(2)-SigP(3)) < tol  && abs(SigP(1)-SigP(3)) < tol)
        
        normal(1,1) = 1.0;
        normal(2,2) = 1.0;
        normal(3,3) = 1.0;
        
        % !	-- two equal eigenvalues -----
    else
        if (abs(SigP(1) - SigP(2)) <= tol)
            hj1 = 1 ; hj2 = 2 ; nr = 3;
        elseif (abs(SigP(1) - SigP(3)) <= tol)
            hj1 = 1 ; nr = 2 ; hj2 = 3;
        else
            nr = 1 ; hj1 = 2 ; hj2 = 3;
        end
        
        % !       --- First principal direction -------------------------
        nl(1) = (Sigma(2)-SigP(nr)) * (Sigma(3)-SigP(nr))- Sigma(6)^2;
        nl(2) = -Sigma(4) * (Sigma(3)-SigP(nr))+ Sigma(6) * Sigma(5);
        nl(3) = Sigma(4) * Sigma(6)- (Sigma(2)-SigP(nr)) * Sigma(5);
        leng = sqrt(nl(1)^2 + nl(2)^2 + nl(3)^2);
        normal(:,nr) = nl/leng ;
        
        % !       --- Second principal direction ------------------------
        if (normal(3,nr) > 0.0  || normal(3,nr) < 0.0)
            nl(1) = 1.0 ; nl(2) = 1.0 ;
            nl(3) = -(normal(1,nr)+normal(2,nr))/normal(3,nr) ;
        elseif (normal(2,nr) > 0.0  || normal(2,nr) < 0.0)
            nl(1) = 1.0 ; nl(3) = 1.0 ;
            nl(2) = -(normal(1,nr)+normal(3,nr))/normal(2,nr) ;
        elseif (normal(1,nr) > 0.0  ||  normal(1,nr) < 0.0)
            nl(2) = 1.0 ; nl(3) = 1.0;
            nl(1) = -(normal(2,nr)+normal(3,nr))/normal(1,nr) ;
        end
        leng = sqrt(nl(1)^2 + nl(2)^2 + nl(3)^2);
        normal(:,hj1) = nl/leng ;
        % !       --- Third principal direction -------------------------
        normal(:,hj2)=cross(normal(:,nr),normal(:,hj1));
        % !			Normal(:,hj2) = cross(normal(:,nr),Normal(:,hj1))
    end
    psi = Normal;
end  %! nsigma == ....
% !----------------------------------------------------------------------------------
end %subroutine PrinDirect !--------------------------------------------------
% !----------------------------------------------------------------------------------
function C=Cross(A,B)
% !--------------------------------------------------------------------------
% ! Calculates the cross product (vector product) between the two
% ! three dimensional vectors A and B
% !
% ! INPUT
% !     A, B  (real) Vectors with three components
% !
% ! OUTPUT
% !     C (real) A vector perpendicular to both A and B
% !
% !------------------------------------------------------------------------
% ! Johan Clausen
% ! Esbjerg Department of Engineering
% ! Aalborg University
% ! July 2005
% !------------------------------------------------------------------------
%
%       real(8), intent(in)		::	A(3), B(3)
%       real(8), intent(out)	::  C(3)
%
% !---------------------------------------------

C(1) = A(2)*B(3) - A(3)*B(2);
C(2) = A(3)*B(1) - A(1)*B(3);
C(3) = A(1)*B(2) - A(2)*B(1);

% !--------------------------------------------------------------------------------
end %subroutine Cross !------------------------------------------------
% !--------------------------------------------------------------------------------

function A=TransMatrix(psi,nsigma,ouplP)
% !--------------------------------------------------------------------
% ! Calculates transformation matrix, A, for stress transformation
% ! based on either an angle (plane and axisymmetry) or direction
% ! cosines (general 3D). The ordering of the shear stresses
% ! must be specified, as it influences how A is built
% !
% ! INPUT
% !  - Name- -type,size -   -- Description --
% !     psi   (real,ar(3,3))  Matrix of direction cosines in general
% !                   three-dimensional stress state or the rotation
% !                   angle otherwise. In 3D we have
% !                   sig_tranformed = psi'*sig*psi, where sig
% !                   is a stress tensor
% !     nsigma  (integer,sc)  Number of stress components, i.e. size of A
% !     ouplP (integer,sc)  Only used in plan situations in which it is similar to
% !                   hoop but for the principal stress vector SigP, i.e.
% !                   it is the position of the out-of-plane principal stress.
% !
% ! OUTPUT
% !     A Transformation matrix, with elements ordered according to hoop and ouplP
% !         A stress is transformed according to ("'" signifies matrix transpose)
% !         Sigma_transformed = inv(A')*Sigma and conversely
% !         Sigma = A'*Sigma_transformed
% !         A strain vector, Epsilon, is tranformed according to
% !         Epsilon_transformed = A*Epsilon and conversely
% !         Epsilon = inv(A)*Epsilon_transformed
% !         See e.g. [Cook, Malkus and Plesha, 1989] for further details
% !--------------------------------------------------------------------
% ! Johan Clausen
% ! Esbjerg Department of Engineering
% ! Aalborg University
% ! July 2005
% !--------------------------------------------------------------------
%
% !-----------------
% integer(4), intent(in) :: nsigma, ouplP
% real(8), intent(in)    :: psi(3,3)
%
% real(8), intent(out)   :: A(nsigma,nsigma)
% !-----------------
%
% real(8) sin_psi, cos_psi, sin_psi2, cos_psi2, sin_2psi
% real(8), dimension(3,3) :: A1, A2, A3, A4, Asmall
% real(8), dimension(3,4) :: Aint
% integer(4), dimension(3):: Hj1, Hj2
% integer(4) Ext(3,4), ExtP(4,3)
% integer(4) i, j
%
% !---- Plane situations --------------------------------------------
% ! Assumes that the stress components are ordered according to
% ! Sigma = [sig_x sig_y sig_z tau_xy]
% !------------------------------------------------------------------
if (nsigma == 4)
    sin_psi  = sin(psi(1,1)) ; cos_psi = cos(psi(1,1));
    sin_psi2 = sin_psi^2     ; cos_psi2 = cos_psi^2;
    sin_2psi = sin(2.0*psi(1,1));
    
    
    Asmall(1,1) = cos_psi2  ; Asmall(1,2) = sin_psi2;
    Asmall(2,1) = sin_psi2  ; Asmall(2,2) = cos_psi2;
    Asmall(3,1) = -sin_2psi ; Asmall(3,2) = sin_2psi;
    
    Asmall(1,3) =  0.5*sin_2psi;
    Asmall(2,3) = -0.5*sin_2psi;
    Asmall(3,3) = cos_psi2 - sin_psi2;
    
    ExtP = 0.0;
    if (ouplP == 2)
        ExtP(1,1) = 1  ;  ExtP(3,2) = 1  ;  ExtP(4,3) = 1;
    elseif (ouplP == 1)
        ExtP(2,1) = 1  ;  ExtP(3,2) = 1  ;  ExtP(4,3) = 1;
    elseif (ouplP == 3)
        ExtP(1,1) = 1  ;  ExtP(2,2) = 1  ;  ExtP(4,3) = 1;
    end
    
    %         ! IMPORTANT! This order of
    Ext = 0.0; %! Ext assumes that the out-of-plane stress i in position 3!!
    Ext(1,1) = 1 ; Ext(2,2) = 1 ; Ext(3,4) = 1;
    
    Aint = Asmall*Ext;
    A = ExtP*Aint;     %! Transformation matrix
    A(ouplP,3) = 1.0;
    
    %         !---- General three dimensional stress state ------------------------
    %         ! Assumes that the stress components are ordered according to
    %         ! Sigma = [sig_x sig_y sig_z tau_xy tau_xz tau_yz]
    %         !------------------------------------------------------------------
elseif (nsigma == 6)
    
    Hj1 = ([1, 3, 2]);
    Hj2 = ([2, 1, 3]);
    
    A1 = psi^2;
    
    for i = 1:3
        for j = 1:3
            A2(i,j) = psi(i,Hj1(j))*psi(i,Hj2(j));
            A3(i,j) = psi(Hj1(i),j)*psi(Hj2(i),j);
            A4(i,j) = psi(Hj1(i),Hj1(j))*psi(Hj2(i),Hj2(j))+ psi(Hj2(i),Hj1(j))*psi(Hj1(i),Hj2(j));
        end
    end
    A2 = 2*A2;
    
    A(1:3,1:3) = A1 ; A(1:3,4:6) = A2;
    A(4:6,1:3) = A3 ; A(4:6,4:6) = A4;
    A = A';
end ; %! nsigma == ....

%     !----------------------------------------------------------------------------------
end %subroutine TransMatrix !-------------------------------------------------
% !----------------------------------------------------------------------------------

function [D,Dinv]=DlinElas(E,nu,nsigma)
% !--------------------------------------------------------------------------
% ! Calculates linear elastic stiffness matrix D. If nsigma = 3, i.e.
% ! only three stress components plane stress is assumed
% !
% ! INPUT
% !     E   (real) Youngs module, module of elasticity
% !     nu    (real) Poisson's ratio
% !     nsigma  (int)  Number of stress components.
% !             nsigma = 4 in plane problems in general
% !             if nsigma = 3 plane stress is assumed
% !             nsigma = 6 in general stress states
% ! OUTPUT
% !     D   (real) Linear elastic stiffness matrix D
% !           Size is dependent on stress type
% !     Dinv  (real) Inverse of D, optional
% !------------------------------------------------------------------------
% ! Johan Clausen
% ! Esbjerg Department of Engineering
% ! Aalborg University
% ! July 2005
% !------------------------------------------------------------------------

% implicit none
%
% integer(4), intent(in) :: nsigma
% real(8), intent(in) :: E, nu
%
% real(8), intent(out) :: D(nsigma,nsigma)
% real(8), intent(out), optional :: Dinv(nsigma,nsigma)
% !-------------------------
% real(8) c, g

D = zeros(nsigma,nsigma);
Dinv = zeros(nsigma,nsigma);
% !----- Plane stress with three stress components -----------------------
if (nsigma == 3)
    
    D(1,1) = 1.0 ; D(1,2) = nu;
    D(2,1) = nu ; D(2,2) = 1.0;
    D(3,3) = (1-nu)/2;
    D = E/(1-nu*nu) * D;
    
    
    Dinv(1,1) = 1.0 ; Dinv(1,2) = -nu;
    Dinv(2,1) = -nu   ; Dinv(2,2) = 1.0;
    Dinv(3,3) = 2.0*(1+nu);
    Dinv = Dinv/E;
    
    %         !----- Plane stress, strain or axisymmetry (four stress components) ---
elseif(nsigma == 4)
    D(1,1) = 1.0 - nu ; D(1,2) = nu         ; D(1,3) = nu;
    D(2,1) = nu         ; D(2,2) = 1.0 - nu ; D(2,3) = nu;
    D(3,1) = nu         ; D(3,2) = nu         ; D(3,3) = 1.0- nu;
    D(4,4) = (1-2.0*nu)/2;
    D = E/((1+nu)*(1-2.0*nu)) * D;
    
    Dinv(1,1) = 1.0 ; Dinv(1,2) = -nu     ; Dinv(1,3) = -nu;
    Dinv(2,1) = -nu   ; Dinv(2,2) = 1.0   ; Dinv(2,3) = -nu;
    Dinv(3,1) = -nu   ; Dinv(3,2) = -nu     ; Dinv(3,3) = 1.0;
    Dinv(4,4) = 2.0*(1+nu);
    Dinv = Dinv/E;
    %             !----- General stress state ---------------------------------------------
elseif(nsigma == 6)
    c = E/((1+nu)*(1-2.0*nu));
    g = E/(2.0*(1+nu));
    
    D(1,1) = (1-nu)*c ; D(1,2) = nu*c   ; D(1,3) = D(1,2);
    D(2,1) = D(1,2)   ; D(2,2) = D(1,1) ; D(2,3) = D(1,2);
    D(3,1) = D(1,2)   ; D(3,2) = D(1,2) ; D(3,3) = D(1,1);
    D(4,4) = g;
    D(5,5) = g;
    D(6,6) = g;
    
    Dinv(1,1) = 1.0 ; Dinv(1,2) = -nu   ; Dinv(1,3) = -nu;
    Dinv(2,1) = -nu   ; Dinv(2,2) = 1.0 ; Dinv(2,3) = -nu;
    Dinv(3,1) = -nu   ; Dinv(3,2) = -nu   ; Dinv(3,3) = 1.0;
    Dinv(4,4) = 2.0*(1+nu);
    Dinv(5,5) = 2.0*(1+nu);
    Dinv(6,6) = 2.0*(1+nu);
    Dinv = Dinv/E;
    
end


end %subroutine  DlinElas !-------------------------------------------------
%             !--------------------------------------------------------------------------------
%%

function [S,I1,J2,J3,lode,sin3lode]=Invariants(Sigma,nsigma)
% !--------------------------------------------------------------------------
% ! Calculates the deviator stress vector of a stress vector Sigma, along
% ! with several stress invariants.
% !
% ! INPUT
% !  - Name -	-- description --
% !	Sigma	(real,ar(nsigma): Stress vector. Size is dependent on the type of stress
% !			  as given by nsigma.
% !			  nsigme = 3: It is assumed that Sigma is principal stresses,
% !						  Sigma = [SigP_1 SigP_2 SigP_3]
% !			  nsigma = 4: Plane situations, i.e. plane stress, plane strain
% !						  and axisymmetry. The ordering of the stresses must
% !						  follow Sigma = [sig_x sig_y sig_z tau_xy] , i.e.
% !						  sig_z is the out-of-plane stress.
% !			  nsigma = 6  General 3D. The stresses must be ordered according
% !						  to Sigma = [sig_x sig_y sig_z tau_xy tau_xz tau_yz]
% !	nsigma	(integer,sc): Number of stress components, i.e. the size of Sigma.
% !
% ! OUTPUT
% !	S		(real, ar(nsigma)): Deviator stress vector. Same size as Sigma
% !	I1		(real,scalar): First stress invariant.
% !	J2		(real,scalar): Second deviator stress invariant.
% !	J3		(real,scalar): Third deviator stress invariant.
% !	lode	(real,scalar): Lode angle in radians. In the
% !			  range -pi/6 <= lode <= pi/6,    (pi/6 = 30 deg)
% !	sin3lode (real,scalar): sin(3*lode)
% !
% !---------------------------------------------------------------------------
% ! Johan Clausen
% ! Deparment of Civil Engineering
% ! Aalborg University
% ! May 2008
% !--------------------------------------------------------------------------
% implicit none
% integer(4), intent(in) :: nsigma
% real(8), intent(in) :: Sigma(nsigma)
%
% real(8), intent(out) :: S(nsigma), I1, J2
% real(8), intent(out) :: J3, lode, sin3lode
%
% !--------------------------------------------------------------------------

I1 = sum(Sigma(1:3));
if (nsigma == 3)  %! Sigma is the principal stresses
    S = Sigma - I1*0.333333333333333333; %! Deviator stress vector
    J2 = 0.5*sum(S^2);
    J3 = S(1)*S(2)*S(3);
    
elseif (nsigma == 4)  %! Plane situation, including axisymmetry
    S(1:3) = Sigma(1:3) - I1*0.333333333333333333; %! Deviator stress vector
    S(4) = Sigma(4);
    J2 = 0.5 * sum(S(1:3)^2) + S(4)^2;
    
    J3 = (S(1)^3 + S(2)^3 + S(3)^3 +3.0*S(4)^2 *(S(1)+S(2)))*0.3333333333333333;
    
elseif (nsigma == 6)  %! General three-dimensional stress state
    S(1:3) = Sigma(1:3) - I1*0.33333333333333333; %! Deviator stress vector
    S(4:6) = Sigma(4:6);
    
    J2 = 0.5 * sum(S(1:3)^2) + sum(S(4:6)^2);
    
    J3 = ( sum(S(1:3)^3) + 6.0*S(4)*S(5)*S(6) +3.0*(S(1)*(S(4)^2 + S(5)^2)...
        +S(2)*(S(4)^2 + S(6)^2) +S(3)*(S(5)^2 + S(6)^2)) )*0.3333333333333333;
    
end

if (J2 > 0.0)
    sin3lode = -2.598076211353316 * J3 / J2^(1.5);
else
    sin3lode = 0.0;
end

if (sin3lode <= -1.0)
    lode = -0.5235987755982988; %! -pi/6
elseif (sin3lode >= 1.0)
    lode = +0.5235987755982988; %! +pi/6
else
    lode = asin(sin3lode)/3.0;
end

%         !----------------------------------------------------------------------------------
end %subroutine Invariants !--------------------------------------------------
%%
function [eps_v, eps_d]=strain_invariant(eps,nsigma) % eps's shear strain components are not multiplied by 2
%
if (nsigma==4)
    eps_v = eps(1) + eps(2) + eps(3);
    e1 = eps(1) - eps_v / 3.0;
    e2 = eps(2) - eps_v / 3.0;
    e3 = eps(3) - eps_v / 3.0;
    e4 = eps(4);
    eps_d =sqrt(2.0/3.0*(e1*e1+e2*e2+e3*e3+2.0*e4*e4));
else
    eps_v = eps(1) + eps(2) + eps(3);
    e1 = eps(1) - eps_v / 3.0;
    e2 = eps(2) - eps_v / 3.0;
    e3 = eps(3) - eps_v / 3.0;
    e4 = eps(4);
    e5 = eps(5);
    e6 = eps(6);
    eps_d =sqrt(2.0/3.0*(e1*e1+e2*e2+e3*e3+2.0*e4*e4+2.0*e5*e5+2.0*e6*e6));
end
end
%%
function [c]=c_reduction(cp,cr,eta,hsv)
c=cr+(cp-cr)*exp(-eta*hsv(18));
end