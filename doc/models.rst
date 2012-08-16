Overview of the models
======================

Each model exists in two forms: a MOMP-only form that can be used to study the
properties of different Bcl-2 reaction topologies as an isolated module, and a
"full TRAIL" form, in which the different MOMP models are embedded in the full
extrinsic apoptosis pathway which begins with TRAIL stimulation.

For the full TRAIL models, the upstream and downstream pathway components
and reaction topologies are re-used from the previously published EARM
1.0. (Albeck et al. 2008 PLoS Biol).

MOMP module "boundaries"
------------------------

In the interest of consistency, all of the models have been defined with the
same boundaries in terms of their position in the overall extrinsic apoptotis
pathway: they are all triggered by the addition of an active BH3-only species
(e.g., tBid) as their most "upstream" event, and they all result in the release
in one or more mitochondrial substances (e.g. Cytochrome C and/or Smac) as
their most downstream event. This represents a compromise between the approach
of the MOMP models described in Albeck et al (in which caspase 8, rather than
tBid, served as the input) and the models of the Shen group, in which active
Bax or Bax pores, rather than Cytochrome C or Smac, served as the output.

While these interface boundaries represent the default condition, they can be
modified by passing parameters in to the module macro. For example, by setting
`do_pore_transport=False` in the call to one of the Shen models, the Cytochrome
C and Smac release reactions are not added, and the models can be directly
compared to their originally published versions. Similarly, the upstream
caspase-8/Bid reactions can be added to the Albeck MOMP models to make them
consistent with their published versions.

Since our purpose in using these models is primarily to embed them in a common
pathway context, rather than to reproduce previous results for posterity, our
conclusion in working with them was that it is better to have a consistent
interface by default and reproduce published results by modifying the model
rather than implement the model as published by default and then have to
specifically modify each one separately to fit the pathway context
appropriately.

The models in EARM
------------------

Here is a list of the models incorporated into EARM. More detailed descriptions
of each model, along with the source code, are found in :doc:`modules/index`.

- EARM 2.0, Embedded
- EARM 2.0, Indirect
- EARM 2.0, Direct
- "Minimal Model" (Figure 11b) from Albeck et al. (2008) [1]_
- "Model B + Bax multimerization" (Figure 11c) from Albeck et al. (2008) [1]_
- "Model C + mitochondrial transport" (Figure 11d) from Albeck et al. (2008)
  [1]_
- "Current model" (Figure 11e) from Albeck et al. (2008) [1]_
- "Current model + cooperativity" (Figure 11f) from Albeck et al. (2008) [1]_
- Deterministic model from Chen et al. (2007), Biophysical Journal [2]_
- Indirect model from Chen et al. (2007), FEBS Letters [3]_
- Direct model from Chen et al. (2007), FEBS Letters [3]_
- Direct model from Cui et al. (2008) [4]_
- Direct model 1 from Cui et al. (2008) [4]_
- Direct model 2 from Cui et al. (2008) [4]_
- Model incorporating Bad phosphorylation from Howells et al. (2011) [5]_

.. [1] Albeck, J. G., Burke, J. M., Spencer, S. L., Lauffenburger, D. A., and
   Sorger, P. K. (2008). Modeling a snap-action, variable-delay switch
   controlling extrinsic cell death. PLoS Biology, 6(12), 2831-2852.
   :doi:`10.1371/journal.pbio.0060299` :pmid:`19053173`.

.. [2] Chen, C., Cui, J., Lu, H., Wang, R., Zhang, S., & Shen,
   P. (2007). Modeling of the role of a Bax-activation switch in the
   mitochondrial apoptosis decision. Biophysical Journal, 92(12),
   4304-4315. :doi:`10.1529/biophysj.106.099606` :pmid:`17400705`.

.. [3] Chen, C., Cui, J., Zhang, W., & Shen, P. (2007). Robustness analysis
   identifies the plausible model of the Bcl-2 apoptotic switch. FEBS letters,
   581(26), 5143-5150. :doi:`10.1016/j.febslet.2007.09.063` :pmid:`17936275`.

.. [4] Cui, J., Chen, C., Lu, H., Sun, T., & Shen, P. (2008). Two independent
   positive feedbacks and bistability in the Bcl-2 apoptotic switch. PLoS ONE,
   3(1), e1469. :doi:`10.1371/journal.pone.0001469` :pmid:`18213378`.

.. [5] Howells, C. C., Baumann, W. T., Samuels, D. C., & Finkielstein,
   C. V. (2010). The Bcl-2-associated death promoter (BAD) lowers the threshold
   at which the Bcl-2-interacting domain death agonist (BID) triggers
   mitochondria disintegration. Journal of Theoretical
   Biology. :doi:`10.1016/j.jtbi.2010.11.040` :pmid:`21130780`.

How to use the models
---------------------

Structure of the python packages. earm, earm.mito. To use a model, run the
Python statement::

    from earm.mito.albeck_11b import model

That's it. You now have a model object that you can query, simulate, perform
parameter estimation on, etc. If you wanted the full TRAIL version, you would
simply run::

    from earm.trail.albeck_11b import model

If you want to work with multiple models at the same time (e.g., to compare
them), you can write::

    from earm.trail.chen2007_indirect import model as indirect
    from earm.trail.chen2007_direct import model as direct

For more information on the kinds of analysis you can do using PySB models,
see the PySB documentation. <TODO>

Parameter values
----------------

Parameter values (both rate constants and initial protein concentrations) are
embedded directly in the model code rather than in a separate table or file.
The values in the model definition represent estimates or nominal values and
can be easily overridden using values obtained (for example) by measurement or
parameter estimation algorithms.  We do not maintain a separate list or table of
parameter values, as we have found that the clearest description of the
meaning of a rate parameter is the macro or rule statement in which it is
embedded.

If desired, lists of all model parameters can be obtained via the parameters
instance variable of the model object, i.e.::

    model.parameters

A list of all parameter names can be obtained using the list comprehension::

    [p.name for p in model.parameters]

The code is meant to be read!
-----------------------------

As much as possible, we have attempted to make the code for models themselves
transparent and well-documented. The documentation for each model topology has
been embedded inline in the model source code, and the documentation provided in
the :doc:`modules/index` section of the documentation is drawn directly from
this source.

Moreover, the models have been written using a high-level vocabulary of
frequently re-used macros, with the aim of revealing broad similarities
and differences between models. The models thus consist of statements such as::

    translocate_tBid_Bax_BclXL()
    catalyze(Bid(state='T'), Bax(state='M'), Bax(state='A'), klist)

which can be read as saying that "tBid, Bax and BclXL translocate [to the
mitochondrial membrane], and tBid catalyzes Bax from a Mitochondrial (but
inactive) state to an Active state." Understanding the precise mechanisms of
these macros (as expressed in terms of rules and reactions) takes some
familiarity with their implementation, but as there is a fairly limited set
of macros, this should hopefully not present a significant barrier.
