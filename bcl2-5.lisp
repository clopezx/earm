(in-package :b-user)
(format t "READING RULES ~%")
(include b-user/ode-biochem)
;;;;
;;;; ------BIOLOGY HERE:
;;;;
;;;;
;;;  A stereotyped cell (used in Birgit's model) which can
;;;  be used to define internalization and other complex reactions
;;;  in which species move between sublocations of the cell.
;;;

;; RATES FILE
(include @FOLDER/bcl2-5_rates)
(format t "DONE READING RATES ~%")
;; cell definition
;; all lengths in um, all areas in um^2 and all volumes in um^3
(defcon cell-location (location)
  (&optional (id := *name*) 
	     &property
	     (cytoplasm compartment :#= [[compartment] {.size.value := 1.0}]); 1000.0, um^3 cell volume of 1pL (in um^3)
	     (dish      compartment :#= [[compartment] {.size.value := 1.0}]); 1.0E6, um^3 1E-9L = dish volume 
	     (mitochondria compartment :#= [compartment])
	     (cell-membrane membrane :#= (let ((d .dish)
					       (c .cytoplasm))
					   [[membrane] :outer d :inner c {.size.value := 1.0}])); 706 um^2, 7.5um cell radius
	     (mito-membrane membrane :#= (let ((e .mitochondria)
					       (c .cytoplasm))
					   [[membrane] :inner e :outer c {.size.value := 1.0}])); 900 um^2, 300 mito, 3um^2/mito
	     (inverse-mito-membrane membrane :#= .mito-membrane.inverse))) ; invert inner and outer
;;
(format t "DONE SETTING CELL ~%")

;; SPECIES DEFINITIONS
;; RECEPTOR LAYER
(defmonomer TRAIL "trail ligand"
  CDREC )
(defmonomer (DRECEPTOR membrane) "Death receptor for TRAIL"
  CTRAIL CCSP8 (asite :states (member :active :inactive) :default :inactive))
(defmonomer caspase8 "caspase 8 with on/off site"
  CSUB CMCL1 CCSP6 CDREC (asite :states (member :active :inactive) :default :inactive))
(defmonomer BIDC "BID, 6 binding sites but only one can bind at a time"
  CCSP8 CPRO CANTI (asite :states (member :active :inactive) :default :inactive)) ; inactive = untruncated
(defmonomer (BIDM membrane) "BID, 6 binding sites but only one can bind at a time"
  CCSP8 CPRO CANTI (asite :states (member :active :inactive) :default :inactive)) ; inactive = untruncated
(defmonomer BAXC "Bax documentation"
  CACTIVATOR CANTI (asite :states (member :active :inactive) :default :inactive))
(defmonomer (BAXM membrane) "Bax documentation"
  CDIM CTETRA CACTIVATOR CANTI (asite :states (member :active :inactive) :default :inactive))
(defmonomer BIM "Bim documentation"
  CPRO CANTI)
(defmonomer (BAK membrane) "Bak documentation"
  CDIM CTETRA CACTIVATOR CANTI (asite :states (member :active :inactive) :default :inactive))
(defmonomer BCLXLC "BCL-XL documentation"
  CBH3 CPRO)
(defmonomer (BCLXLM membrane) "BCL-XL documentation"
  CBH3 CPRO)
(defmonomer (BCL2 membrane) "BCL2 documentation"
  CBH3 CPRO)
;; Added state to MCL1 to account for degradation. Default of MCL1 is active (unless specified).
(defmonomer MCL1C "MCL1 documentation" 
  CCSP8 CBH3 CPRO (asite :states (member :active :degraded) :default :active))
(defmonomer (MCL1M membrane) "MCL1 documentation" 
  CCSP8 CBH3 CPRO (asite :states (member :active :degraded) :default :active))
(defmonomer BAD "BAD documentation"
  CANTI)
(defmonomer NOXA "NOXA documentation"
  CANTI)
(defmonomer MULE "MULE documentation"
  CANTI)
(defmonomer caspase3 "caspase 3"
  CCSP8 CCSP6 CMCL1 CXIAP CPARP CAPOP
  (asite :states (member :active :inactive) :default :inactive)  ;activation site
  (dsite :states (member :ubiq   :nonubiq ) :default :nonubiq )) ;ubiquination site
(defmonomer caspase6 "caspase 6"
  CCSP3 CCSP8 (asite :states (member :active :inactive) :default :inactive))
(defmonomer CYTC "cytochrome C documentation"
  CPORE CAPAF)
(defmonomer SMAC "SMAC/DIABLO documentation"
  CPORE CXIAP)
(defmonomer (BAXMPORE membrane) "pore"
  CMITOP)
(defmonomer (BAKPORE  membrane) "pore"
  CMITOP)
(defmonomer XIAP "x-linked inhibitor of apoptosis proteins"
  CCSP3 CAPOP CSMAC)
(defmonomer PARP "Poly(ADP-ribose) polymerase"
  CCSP3 (asite :states (member :active :inactive) :default :active))
(defmonomer APAF "Apoptosis activating factor 1"
  CCYTC CCSP9 (asite :states (member :active :inactive) :default :inactive))
(defmonomer caspase9 "caspase 9"
  CAPAF)
(defmonomer apoptosome "Apoptosome"
  CCSP3 CXIAP)

(format t "DONE SETTING MONOMERS ~%")

;;; ************************************************************************
;;; RULEZ
;;; ************************************************************************

;;; TRAIL BINDING TO RECEPTOR
[{[TRAIL CDREC._ __] @ :outer + [DRECEPTOR CTRAIL._ __] <<->> [[TRAIL CDREC.1][DRECEPTOR CTRAIL.1]]}
 (.set-rate-function 'mass-action :fwd KTRAILRECF :rev KTRAILRECR)]
[{[[TRAIL CDREC.1][DRECEPTOR CTRAIL.1]] ->> [DRECEPTOR CTRAIL._ asite.active]} ;NOTICE TRAIL DISAPPEARS
 (.set-rate-function 'mass-action KTRAILRECC)]

;;; EXCLUDING FLIP INHIBITOR OF DEATH RECEPTOR

;;; CSP 8 activation
[{[DRECEPTOR CCSP8._ asite.active] + [caspase8 CDREC._ asite.inactive] @ :inner <<->>
  [[DRECEPTOR CCSP8.1 asite.active][caspase8 CDREC.1 asite.inactive]]}
 (.set-rate-function 'mass-action :fwd KDRECCSP8F :rev KDRECCSP8R)]
[{[[DRECEPTOR CCSP8.1 asite.active][caspase8 CDREC.1 asite.inactive]] ->>
  [DRECEPTOR CCSP8._ asite.active] + [caspase8 CDREC._ asite.active] @ :inner }
 (.set-rate-function 'mass-action KDRECCSP8C)]

;;; MCL1 degradation
;; MCL1 gets degraded in the presence of cycloheximide
;[{[MCL1C asite.active] ->> [MCL1C asite.degraded]}
;  (.set-rate-function 'mass-action KMCL1CAMCL1CD)]
;[{[MCL1M asite.active] ->> [MCL1M asite.degraded]}
;  (.set-rate-function 'mass-action KMCL1MAMCL1MD)]

;;; BH3 ACTIVATORS BINDING BAX BAK
;; \wikidoc{bh3_activators_bax_bak}
;; CYTOSOL 
;; can happen both in the cytoplasm and in the membrane
(with-data-table 
    (:rows ($R1 @A) :cols $R2 :cells ($Kf $Kr $Kc @B @C) :ignore _)
  ((                       BAXC                                              BAK                )
   ((BIDC (asite.active)) (KBIDCBAXCF KBIDCBAXCR KBIDCBAXCC NIL NIL)        (KBIDCBAKF KBIDCBAKR KBIDCBAKC (@ :outer) NIL))
   ((BIDM (asite.active)) (KBIDMBAXCF KBIDMBAXCR KBIDMBAXCC NIL (@ :outer)) (KBIDMBAKF KBIDMBAKR KBIDMBAKC  NIL NIL))
   ((BIM   NIL          ) (KBIMBAXCF  KBIMBAXCR  KBIMBAXCC  NIL NIL)        (KBIMBAKF  KBIMBAKR  KBIMBAKC  (@ :outer) NIL)))
  ;; FORWARD (BINDING)
  [{[$R1 CPRO._ @A] @B + [$R2 CACTIVATOR._ asite.inactive] @C  <<->> [[$R1 CPRO.1 @A][$R2 CACTIVATOR.1 asite.inactive]]}
   (.set-rate-function 'mass-action :fwd $Kf  :rev $Kr)]
  [{[[$R1 CPRO.1 @A][$R2 CACTIVATOR.1 asite.inactive]] ->> [$R1 CPRO._ @A] @B + [$R2 CACTIVATOR._ asite.active] @C}
   (.set-rate-function 'mass-action $Kc)])

;; Translocations
;;\wikidoc{translocation}
(with-substitution-table 
 (($R1     @A             $R2     $Kf            $Kr  ) 
  (BIDC   (asite.active)  BIDM    KBIDCBIDMF     KBIDCBIDMR )
  (BAXC   (asite.active)  BAXM    KBAXCBAXMF     KBAXCBAXMR )
  (BCLXLC  NIL            BCLXLM  KBCLXLCBCLXLMF KBCLXLCBCLXLMR)
  (MCL1C   NIL            MCL1M   KMCL1CMCL1MF   KMCL1CMCL1MR  ))
 [{[$R1 @A] @ cell-location.cytoplasm <<->> [$R2 @A] @ :mito-membrane}
  (.set-rate-function 'mass-action :fwd $Kf :rev $Kr)])

;; AUTO-ACTIVATION
;; BAK/BAX "spontaneously" activates itself (australian model)
(with-substitution-table
 (($R1   $Kf            $Kr                         )
  (BAXC  KBAXAUTOACTF   KBAXAUTOACTR                )
  (BAK   KBAKAUTOACTF   KBAKAUTOACTR                ))
  [{[$R1 asite.inactive __] <<->> [$R1 asite.active __]}
 (.set-rate-function 'mass-action :fwd  $Kf :rev $Kr)])

;; BAK/BAX Oligomerization (dimer formation)
(with-substitution-table
 (($R1     $Kf          $Kr  )
  (BAXM    KBAXMDIMF   KBAXMDIMR)
  (BAK     KBAKDIMF    KBAKDIMR ))
  [{[$R1 CDIM._ asite.active __] + [$R1 CDIM._ asite.active __]  <<->> 
    [[$R1 CDIM.1 asite.active __][$R1 CDIM.1 asite.active __]]}
   (.set-rate-function 'mass-action :fwd $Kf :rev $Kr)])

;; BAK/BAXM Oligomerization (tetramer formation)
;; littleb gets freaked out if you allow multibonds in this situation
;; through a tetramer. I replaced the tetramer with 
(with-substitution-table
 (($R1   $Kf         $Kr )
  (BAXM  KBAXMTETF  KBAXMTETR)
  (BAK   KBAKTETF   KBAKTETR ))
  [{[[$R1 CDIM.1 asite.active __][$R1 CDIM.1 asite.active __]] + 
    [[$R1 CDIM.2 asite.active __][$R1 CDIM.2 asite.active __]]  <<->>
    [[$R1 CDIM.1 CTETRA.3 asite.active __][$R1 CDIM.1 CTETRA.4 asite.active __]
     [$R1 CDIM.2 CTETRA.3 asite.active __][$R1 CDIM.2 CTETRA.4 asite.active __]]}
   (.set-rate-function 'mass-action :fwd $Kf :rev $Kr )])

;; BAK/BAXM tetramer becomes a "pore"
(with-substitution-table
 (($R1   $R2        $Kf )
  (BAXM  BAXMPORE   KBAXMPORE)
  (BAK   BAKPORE    KBAKPORE ))
  [{[[$R1 CDIM.1 CTETRA.3 asite.active __][$R1 CDIM.1 CTETRA.4 asite.active __]
     [$R1 CDIM.2 CTETRA.3 asite.active __][$R1 CDIM.2 CTETRA.4 asite.active __]] 
    ->> [$R2 __] }
   (.set-rate-function 'mass-action $Kf )])

;; ANTIAPOPTOTIC inhibiting BAX, BAK
(with-data-table (:rows ($R1 @A) :cols ($R2 @B) :cells ($Kf $Kr @C @D) :ignore _)
  ((                      (BAXC NIL)    (BAXM (CDIM.* CTETRA.*))     (BAK (CDIM.* CTETRA.*)))
   ((BCLXLC NIL)          (KBCLXLCBAXCF KBCLXLCBAXCR NIL NIL) 
	                                (KBCLXLCBAXMF KBCLXLCBAXMR (@ :outer) NIL) 
	                                                    (KBCLXLCBAKF  KBCLXLCBAKR (@ :outer) NIL))
   ((BCLXLM NIL)          (KBCLXLMBAXCF KBCLXLMBAXCR NIL (@ :outer))      
                                        (KBCLXLMBAXMF KBCLXLMBAXMF NIL NIL)      
                                                            (KBCLXLMBAKF KBCLXLMBAKR NIL NIL))
   ((BCL2   NIL)          (KBCL2BAXCF KBCL2BAXCR NIL (@ :outer))      
                                        (KBCL2BAXMF KBCL2BAXMF NIL NIL)      
                                                            (KBCL2BAKF KBCL2BAKR NIL NIL))
   ((MCL1C  (asite.active)) (KMCL1CBAXCF KMCL1CBAXCR NIL NIL)      
                                        (KMCL1CBAXMF KMCL1CBAXMF (@ :outer) NIL)      
                                                            (KMCL1CBAKF KMCL1CBAKR (@ :outer) NIL))
   ((MCL1M  (asite.active)) (KMCL1MBAXCF KMCL1MBAXCR NIL (@ :outer))      
                                        (KMCL1MBAXMF KMCL1MBAXMF NIL NIL)      
                                                            (KMCL1MBAKF KMCL1MBAKR NIL NIL)))
  [{[$R1 @A CPRO._ __] @C + [$R2 @B CANTI._ asite.*  __] @D <<->> [[$R1 @A CPRO.1 __][$R2 @B CANTI.1 asite.*  __]]}
   (.set-rate-function 'mass-action :fwd $Kf :rev $Kr)])
  
;; ANTIAPOPTOTIC binding BH3s
(with-data-table (:rows ($R1 @A) :cols ($R2 @B) :cells ($Kf $Kr @C @D) :ignore _)
  ((                       (BIDC (asite.active)) (BIDM (asite.active)) (BIM NIL) (BAD NIL) (NOXA NIL) (MULE NIL))
   ((BCLXLC NIL)           (KBCLXLCBIDCF KBCLXLCBIDCR NIL        NIL)  ; with BIDC
                           (KBCLXLCBIDMF KBCLXLCBIDMR (@ :outer) NIL)  ; with BIDM
                           (KBCLXLCBIMF  KBCLXLCBIMR  NIL        NIL)  ; with BIM
                           (KBCLXLCBADF  KBCLXLCBADR  NIL        NIL)  ; with BAD
                            _                                          ; with NOXA
                            _                                        ) ; with MULE
   ((BCLXLM NIL)           (KBCLXLMBIDCF KBCLXLMBIDCR NIL (@ :outer))  ; with BIDC
                           (KBCLXLMBIDMF KBCLXLMBIDMR NIL        NIL)  ; with BIDM
                           (KBCLXLMBIMF  KBCLXLMBIMR  NIL (@ :outer))  ; with BIM
                           (KBCLXLMBADF  KBCLXLMBADR  NIL (@ :outer))  ; with BAD
                            _                                          ; with NOXA
                            _                                        ) ; with MULE
   ((BCL2   NIL)           (KBCL2BIDCF   KBCL2BIDCR   NIL (@ :outer))  ; with BIDC
                           (KBCL2BIDMF   KBCL2BIDMR   NIL        NIL)  ; with BIDM
                           (KBCL2BIMF    KBCL2BIMR    NIL (@ :outer))  ; with BIM
                           (KBCL2BADF    KBCL2BADR    NIL (@ :outer))  ; with BAD
                            _                                          ; with NOXA
                            _                                        ) ; with MULE
   ((MCL1C (asite.active)) (KMCL1CBIDCF  KMCL1CBIDCR  NIL        NIL)  ; with BIDC
                           (KMCL1CBIDMF  KMCL1CBIDMR  (@ :outer) NIL)  ; with BIDM
                           (KMCL1CBIMF   KMCL1CBIMR   NIL        NIL)  ; with BIM
                            _                                          ; with BAD
			   (KMCL1CNOXAF  KMCL1CNOXAR  NIL        NIL)  ; with NOXA
                           (KMCL1CMULEF  KMCL1CMULER  NIL        NIL)) ; with MULE
   ((MCL1M (asite.active)) (KMCL1MBIDCF  KMCL1MBIDCR  NIL (@ :outer))  ; with BIDC
                           (KMCL1MBIDMF  KMCL1MBIDMR  NIL        NIL)  ; with BIDM
                           (KMCL1MBIMF   KMCL1MBIMR   NIL (@ :outer))  ; with BIM
                            _                                          ; with BAD
			   (KMCL1MNOXAF  KMCL1MNOXAR  NIL (@ :outer))  ; with NOXA
                           (KMCL1MMULEF  KMCL1MMULER  NIL (@ :outer)))); with MULE
  [{[$R1 @A CBH3._ __] @C + [$R2 @B CANTI._ __] @D <<->> [[$R1 @A CBH3.1 __][$R2 CANTI.1 @B __]]}
   (.set-rate-function 'mass-action :fwd $Kf :rev $Kr)])

;; MCL1 degradation  
(with-data-table (:rows $R1 :cols $R2 :cells ($Kc @A @B) :ignore _)
  ((         NOXA                        MULE)
   (MCL1C   (KMCL1CNOXAC NIL  NIL)        (KMCL1CMULEC NIL  NIL))
   (MCL1M   (KMCL1MNOXAC NIL (@ :outer))  (KMCL1MMULEC NIL (@ :outer))))
  [{[[$R1 asite.active CBH3.1 __][$R2 CANTI.1 __]] ->> [$R1 asite.degraded CBH3._ __] @A + [$R2 CANTI._ __] @B }
   (.set-rate-function 'mass-action $Kc)])


;; CSP8 activates BIDC and CSP3
(with-data-table (:rows $R1 :cols $R2 :cells ($Kf $Kr $Kc) :ignore _)
  ((              caspase8                                          )
   (BIDC         (KBIDCCSP8F KBIDCCSP8R KBIDCCSP8C  ))
   (caspase3     (KCSP3CSP8F KCSP3CSP8R KCSP3CSP8C)))
  [{[$R1 asite.inactive CCSP8._  __] + [$R2 asite.active CSUB._ __]  <<->> [[$R1 asite.inactive CCSP8.1][$R2 asite.active CSUB.1]]}
   (.set-rate-function 'mass-action :fwd $Kf  :rev $Kr)]
  [{[[$R1 asite.inactive CCSP8.1][$R2 asite.active CSUB.1]] ->> [$R1 asite.active CCSP8._ ] + [$R2 asite.active CSUB._]}
   (.set-rate-function 'mass-action $Kc)])

;; CSP3 activates CSP6
[{[caspase3 CCSP6._ asite.active] + [caspase6 CCSP3._ asite.inactive] <<->> 
  [[caspase3 CCSP6.1 asite.active][caspase6 CCSP3.1 asite.inactive]]}
 (.set-rate-function 'mass-action :fwd KCSP3CSP6F :rev KCSP3CSP6R)]
[{[[caspase3 CCSP6.1 asite.active][caspase6 CCSP3.1 asite.inactive]] ->> 
  [caspase3 CCSP6._ asite.active] + [caspase6 CCSP3._ asite.active]}
 (.set-rate-function 'mass-action KCSP3CSP6C)]

;; CSP3 gets tagged for ubiquination by XIAP
[{[caspase3 CXIAP._ asite.active dsite.nonubiq] + [XIAP CCSP3._ ] <<->> 
  [[caspase3 CXIAP.1 asite.active dsite.nonubiq][XIAP CCSP3.1]]}
 (.set-rate-function 'mass-action :fwd KCSP3XIAPF :rev KCSP3XIAPR)]
[{[[caspase3 CXIAP.1 asite.inactive dsite.nonubiq][XIAP CCSP3.1]] ->> 
  [XIAP CCSP3._] + [caspase3 CXIAP._ asite.inactive dsite.ubiq]}
 (.set-rate-function 'mass-action KCSP3XIAPC)]

;; CSP3 cleaves PARP and inactivates it
[{[caspase3 CPARP._ asite.active] + [PARP CCSP3._] <<->> 
  [[caspase3 CPARP.1 asite.active][PARP CCSP3.1]]}
 (.set-rate-function 'mass-action :fwd KCSP3PARPF :rev KCSP3PARPR)]
[{[[caspase3 CPARP.1 asite.active][PARP CCSP3.1]] ->>
  [caspase3 CPARP._ asite.active] + [PARP CCSP3._ asite.inactive]}
 (.set-rate-function 'mass-action KCSP3PARPC)]

;; CSP6 activates CSP8
[{[caspase6 CCSP8._ asite.active] + [caspase8 CCSP6._ asite.inactive] <<->> 
  [[caspase6 CCSP8.1 asite.active][caspase8 CCSP6.1 asite.inactive]]}
 (.set-rate-function 'mass-action :fwd KCSP6CSP8F :rev KCSP6CSP8R)]
[{[[caspase6 CCSP8.1 asite.active][caspase8 CCSP6.1 asite.inactive]] ->> 
  [caspase6 CCSP8._ asite.active] + [caspase8 CCSP6._ asite.active]}
 (.set-rate-function 'mass-action KCSP6CSP8C)]

;; \wikidoc{CSP8_inhibits_MCL1}
;; caspase8/caspase3 induce  MCL1 degradation
(with-data-table (:rows $R1 :cols $R2 :cells ($Kf $Kr $Kc @A) :ignore _)
  ((           caspase8                           caspase3                        )
   (MCL1M    (KMCL1MCSP8F KMCL1MCSP8R KMCL1MCSP8C (@ :outer)) (KMCL1MCSP3F KMCL1MCSP3R KMCL1MCSP3C (@ :outer)))
   (MCL1C    (KMCL1CCSP8F KMCL1CCSP8R KMCL1CCSP8C  NIL      ) (KMCL1CCSP3F KMCL1CCSP3R KMCL1CCSP3C  NIL)))
    [{[$R1 CCSP8._ asite.active] + [$R2 asite.active CMCL1._] @A  <<->> [[$R1 CCSP8.1 asite.active][$R2 asite.active CMCL1.1]]}
     (.set-rate-function 'mass-action :fwd $Kf :rev $Kr)]
    [{[[$R1 CCSP8.1 asite.active][$R2 asite.active CMCL1.1]] ->> [$R1 CCSP8._ asite.degraded] + [$R2 asite.active CMCL1._] @A}
     (.set-rate-function 'mass-action $Kc)])

;; \wikidoc{BAXP_releases_CYTC}
;; CYTOCHROME C/SMAC released via BAXMPORE/BAKPORE
(with-data-table (:rows $R1 :cols $R2 :cells ($Kf $Kr $Kc) :ignore _)
  ((          BAXMPORE                            BAKPORE )
   (CYTC     (KCYTCBAXPF KCYTCBAXPR KCYTCBAXPC) (KCYTCBAKPF KCYTCBAKPR KCYTCBAKPC))
   (SMAC     (KSMACBAXPF KSMACBAXPR KSMACBAXPC) (KSMACBAKPF KSMACBAKPR KSMACBAKPC)))
  [{[$R1 CPORE._ __] @ :inner + [$R2 CMITOP._ __]    <<->> [[$R1 CPORE.1 __][$R2 CMITOP.1 __]]}
   (.set-rate-function 'mass-action :fwd $Kf :rev $Kr )]
  [{[[$R1 CPORE.1 __][$R2 CMITOP.1 __]]  ->>  [$R1 CPORE._ __] @ :outer  + [$R2 CMITOP._ __]}
   (.set-rate-function 'mass-action $Kc)])

;; CYTOCHROME C complexes with APAF
[{[CYTC CAPAF._] @ cell-location.cytoplasm  + [APAF CCYTC._ asite.inactive] @ cell-location.cytoplasm <<->> 
  [[CYTC CAPAF.1][APAF CCYTC.1 asite.inactive]] @ cell-location.cytoplasm}
 (.set-rate-function 'mass-action :fwd KCYTCAPAFF  :rev KCYTCAPAFR)]

;; APAF and CSP9 form the apoptosome
[{[[CYTC CAPAF.1][APAF CCYTC.1 asite.inactive]] @ cell-location.cytoplasm + [caspase9 CAPAF._] @ cell-location.cytoplasm <<->> [apoptosome] @ cell-location.cytoplasm } 
 (.set-rate-function 'mass-action :fwd KAPAFCSP9F :rev KAPAFCSP9R)]
 
;; apoptosome activates caspase3
[{[apoptosome CCSP3._] + [caspase3 CAPOP._ asite.inactive] <<->>
  [[apoptosome CCSP3.1][caspase3 CAPOP.1 asite.inactive]]}
 (.set-rate-function 'mass-action :fwd KAPOPCSP3F :rev KAPOPCSP3R)]
[{[[apoptosome CCSP3.1][caspase3 CAPOP.1 asite.inactive]] ->>
  [apoptosome CCSP3._] + [caspase3 CAPOP._ asite.active]}
 (.set-rate-function 'mass-action  KAPOPCSP3C)]

;; apoptosome binds to XIAP
[{[apoptosome CXIAP._] + [XIAP CAPOP._] <<->> [[apoptosome CXIAP.1][XIAP CAPOP.1]]}
(.set-rate-function 'mass-action :fwd KAPOPXIAPF :rev KAPOPXIAPR)]

;; smac binds to XIAP
[{[SMAC CXIAP._] + [XIAP CSMAC._] <<->> [[SMAC CXIAP.1][XIAP CSMAC.1]]}
(.set-rate-function 'mass-action :fwd KSMACXIAPF :rev KSMACXIAPR)]


;;----------------------------------------------------------------------------------------
;; INITIAL CONDITIONS
;;----------------------------------------------------------------------------------------

;; DEFINE CELL AND SPECIES IN COMPARTMENTS
(define cell [cell-location])
cell.dish.(contains [TRAIL])
cell.cell-membrane.(contains [DRECEPTOR])
cell.cytoplasm.(contains [caspase8] [caspase3] [caspase6] [caspase9] [BIDC] [BAXC] [BIM] 
			 [BCLXLC]   [MCL1C]    [NOXA]     [BAD]      [MULE] [XIAP] [PARP] 
			 [APAF]     [apoptosome])
cell.mito-membrane.(contains [BIDM] [BAXM] [BAK] [BCLXLM] [BCL2] [MCL1M] [BAXMPORE] [BAKPORE])
cell.mitochondria.(contains [CYTC] [SMAC])

;; INITIAL CONCENTRATIONS/AMOUNTS
{[TRAIL].(in cell.dish).conc.t0 := 3000} ;; molecules/cell original: 3000
{[DRECEPTOR].(in cell.cell-membrane).conc.t0 := 200} ;; molecules/cell
{[caspase8].(in cell.cytoplasm).conc.t0 := 2E4 } ;; molecules/cell
{[caspase3].(in cell.cytoplasm).conc.t0 := 1E4} ;; molecules/cell
{[caspase6].(in cell.cytoplasm).conc.t0 := 1E4 } ;; molecules/cell
{[caspase9].(in cell.cytoplasm).conc.t0 := 1E5 } ;; molecules/cell
{[BIDC].(in cell.cytoplasm).conc.t0 := 4E4 } ;; molecules/cell
{[BIDM].(in cell.mito-membrane).conc.t0 := 0 } ;; molecules/cell
{[BAXC].(in cell.cytoplasm).conc.t0 := 1.0E5 } ;; molecules/cell
{[BAXM].(in cell.mito-membrane).conc.t0 := 0} ;; molecules/cell
{[BIM].(in cell.cytoplasm).conc.t0 := 0} ;; molecules/cell 5E4
{[BAK].(in cell.mito-membrane).conc.t0 := 5E4} ;; molecules/cell
{[BCLXLC].(in cell.cytoplasm).conc.t0 := 1.98E4} ;; molecules/cell
{[BCLXLM].(in cell.mito-membrane).conc.t0 := 200} ;; molecules/cell
{[BCL2].(in cell.mito-membrane).conc.t0 := 2E4} ;; molecules/cell
{[MCL1C].(in cell.cytoplasm).conc.t0 := 1.98E4} ;; molecules/cell
{[MCL1M].(in cell.mito-membrane).conc.t0 := 200} ;; molecules/cell
{[BAD].(in cell.cytoplasm).conc.t0 := 5E4} ;; molecules/cell
{[NOXA].(in cell.cytoplasm).conc.t0 := 5E4} ;; molecules/cell
{[MULE].(in cell.cytoplasm).conc.t0 := 1E4} ;; molecules/cell
{[XIAP].(in cell.cytoplasm).conc.t0 := 1E5} ;; molecules/cell
{[PARP].(in cell.cytoplasm).conc.t0 := 1E6} ;; molecules/cell
{[APAF].(in cell.cytoplasm).conc.t0 := 1E5} ;; molecules/cell
{[apoptosome].(in cell.cytoplasm).conc.t0 := 0} ;; molecules/cell
{[CYTC].(in cell.mitochondria).conc.t0 := 5E5} ;; molecules/cell
{[SMAC].(in cell.mitochondria).conc.t0 := 1E5} ;; molecules/cell
{[BAXMPORE].(in cell.mito-membrane).conc.t0 := 0} ;; molecules/cell
{[BAKPORE].(in cell.mito-membrane).conc.t0 := 0} ;; molecules/cell

;; \wikidoc{create_model}
;(include b/matlab/ode-translation)
;(create-ode-model "input/erb6_rcpt" :ode-comments nil :overwrite t)
(include b/numerica/ode-translation)
(setf *NUMERICA-RATE-STRING-MAX-LENGTH* NIL)
(create-numerica-model "input/bcl2_5" :ode-comments nil :overwrite t :vars (query species.moles).alpha-order
		       :sim-options '("CSVOUTPUT := TRUE") :sim-steps 56600)
;(create-numerica-model "input/bcl2_4" :ode-comments nil :overwrite t :vars (query species.conc)
;		       :sim-options '("CSVOUTPUT := TRUE" "DYNAMIC_REPORTING_INTERVAL := 600.0") :sim-steps 1800)
