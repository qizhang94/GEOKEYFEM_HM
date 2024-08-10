function [FORCE, STIFF] = cntelm2d(OMEGAN, OMEGAT, CFRI, ELXY, ELXYP, LTAN)
%***********************************************************************
% SEARCH CONTACT POINT AND RETURN STIFFNESS AND RESIDUAL FORCE 
% IF CONTACTED FOR NORMAL CONTACT
%***********************************************************************
%
  ZERO = 0.D0; 
  ONE = 1.D0; 
  EPS = 1.E-6;
  P05 = 0.05;
  FORCE=[];
  STIFF=[];
  XT = ELXY(:,3)-ELXY(:,2); 
  XLEN = norm(XT);
  if XLEN < EPS ,return; end
  XTP = ELXYP(:,3)-ELXYP(:,2); XLENP  = norm(XTP);
%
% UNIT NORMAL AND TANGENTIAL VECTOR
  XT = XT/XLEN;
  XTP = XTP/XLENP;
  XN = [-XT(2); XT(1)];
%
% NORMAL GAP FUNCTION Gn = (X_s - X_1).N
  GAPN = (ELXY(:,1)-ELXY(:,2))'*XN;
%
% CHECK IMPENETRATION CONDITION
  if (GAPN >= ZERO) || (GAPN <= -XLEN), return; end
%
% NATURAL COORDINATE AT CONTACT POINT
  ALPHA = (ELXY(:,1) - ELXY(:,2))'*XT/XLEN;
  ALPHA0 = ((ELXYP(:,1)-ELXYP(:,2))'*XTP)/XLENP;
%
% OUT OF SEGMENT
  if (ALPHA > ONE+P05) || (ALPHA < -P05), return; end
%
% CONTACT OCCURS IN THIS SEGMENT
  XLAMBN = -OMEGAN*GAPN;
  XLAMBT = 0;
  LFRIC = 1; if CFRI == 0, LFRIC = 0; end
  if LFRIC
    GAPT = (ALPHA - ALPHA0)*XLENP;
    XLAMBT = -OMEGAT*GAPT;
    FRTOL  = XLAMBN*CFRI;
    LSLIDE = 0;
    if abs(XLAMBT) > FRTOL
      LSLIDE = 1;
      XLAMBT = -FRTOL*sign(GAPT);
    end
  end
%
% DEFINE VECTORS
  NN = [XN; -(ONE-ALPHA)*XN; -ALPHA*XN];
  TT = [XT; -(ONE-ALPHA)*XT; -ALPHA*XT];
  PP = [ZERO; ZERO; -XN; XN];
  QQ = [ZERO; ZERO; -XT; XT];
  CN = NN - GAPN*QQ/XLEN;
  CT = TT + GAPN*PP/XLEN;
%
% CONTACT FORCE
  FORCE = XLAMBN*CN + XLAMBT*CT;
%
% FORM STIFFNESS
  if LTAN
    STIFF = OMEGAN*(CN*CN');
    if LFRIC
      TMP1 = -CFRI*OMEGAN*sign(GAPT);
      TMP2 = -XLAMBT/XLEN;
      if LSLIDE
        STIFF = STIFF + TMP1*(CT*CN') + TMP2*(CN*PP'+PP*CN'-CT*QQ'-QQ*CT');
      else
        STIFF = STIFF + OMEGAT*(CT*CT') + TMP2*(CN*PP'+PP*CN'-CT*QQ'-QQ*CT');
      end
    end
  end
end