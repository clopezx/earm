#-------------------------------------------------------------------------
Model()

# Components, macros, etc.


#-------------------------------------------------------------------------
import base_model

Model(base=base_model.model)

# Add new components
# Modify existing components

#-------------------------------------------------------------------------
def ligand_to_disc_module():
    util.alias_model_components()
    # Components, macros, etc.

def pore_to_parp_module():
    util.alias_model_components()
    # Components, macros, etc.

#-------------------------------------------------------------------------
import earm_modules

Model()

earm_modules.ligand_to_disc_module()

# Implement particular MOMP model
# Components, macros, etc. 

earm_modules.pore_to_parp_module()

#-------------------------------------------------------------------------
class BaseBuilder(object):
    def tBid_activates_Bax(self):
        # Components, macros, etc.

    def Bax_forms_pores(self):
        # Components, macros, etc.

    def build_my_model1(self):
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax()

    def build_my_model2(self):
        self.translocate_tBid_Bax()
        self.tBid_activates_Bax()
        self.Bax_forms_pores()

#-------------------------------------------------------------------------
from base_builder import BaseBuilder

class TwoCptBuilder(BaseBuilder):
    def translocate_tBid_Bax(self):
        # Components, macros, etc.

#-------------------------------------------------------------------------
from base_builder import BaseBuilder

class MultiCptBuilder(BaseBuilder):
    def translocate_tBid_Bax(self):
        # Components, macros, etc.
