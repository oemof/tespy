.. tab-set::

    .. tab-item:: components

        .. rubric:: Component

        Class documentation and example: :py:class:`Component <tespy.components.component.Component>`

        Table of constraints

        =====================  =============  =======================================================================================================================
        Parameter              Description    Method
        =====================  =============  =======================================================================================================================
        mass_flow_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        =====================  =============  =======================================================================================================================



    .. tab-item:: basics

        .. rubric:: CycleCloser

        Class documentation and example: :py:class:`CycleCloser <tespy.components.basics.cycle_closer.CycleCloser>`

        Table of constraints

        ============================  =============  =======================================================================================================================
        Parameter                     Description    Method
        ============================  =============  =======================================================================================================================
        pressure_equality_constraint  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        enthalpy_equality_constraint  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        ============================  =============  =======================================================================================================================

        Table of parameters

        ===============  =============  ============  ============
        Parameter        Description    Quantity      Method
        ===============  =============  ============  ============
        mass_deviation   :code:`None`   mass_flow     :code:`None`
        fluid_deviation  :code:`None`   :code:`None`  :code:`None`
        ===============  =============  ============  ============




        .. rubric:: Sink

        Class documentation and example: :py:class:`Sink <tespy.components.basics.sink.Sink>`




        .. rubric:: Source

        Class documentation and example: :py:class:`Source <tespy.components.basics.source.Source>`

        Table of constraints

        =====================  =============  =======================================================================================================================
        Parameter              Description    Method
        =====================  =============  =======================================================================================================================
        mass_flow_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        =====================  =============  =======================================================================================================================




        .. rubric:: SubsystemInterface

        Class documentation and example: :py:class:`SubsystemInterface <tespy.components.basics.subsystem_interface.SubsystemInterface>`

        Table of constraints

        =====================  =============  =======================================================================================================================
        Parameter              Description    Method
        =====================  =============  =======================================================================================================================
        mass_flow_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        pressure_constraints   :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        enthalpy_constraints   :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        =====================  =============  =======================================================================================================================

        Table of parameters

        ===========  =============  ============  ============
        Parameter    Description    Quantity      Method
        ===========  =============  ============  ============
        num_inter    :code:`None`   :code:`None`  :code:`None`
        ===========  =============  ============  ============



    .. tab-item:: combustion

        .. rubric:: CombustionChamber

        Class documentation and example: :py:class:`CombustionChamber <tespy.components.combustion.base.CombustionChamber>`

        Table of constraints

        ============================  =============  =========================================================================================================================================
        Parameter                     Description    Method
        ============================  =============  =========================================================================================================================================
        mass_flow_constraints         :code:`None`   :py:meth:`mass_flow_func <tespy.components.combustion.base.CombustionChamber.mass_flow_func>`
        reactor_pressure_constraints  :code:`None`   :py:meth:`combustion_pressure_structure_matrix <tespy.components.combustion.base.CombustionChamber.combustion_pressure_structure_matrix>`
        stoichiometry_constraints     :code:`None`   :py:meth:`stoichiometry_func <tespy.components.combustion.base.CombustionChamber.stoichiometry_func>`
        energy_balance_constraints    :code:`None`   :py:meth:`energy_balance_func <tespy.components.combustion.base.CombustionChamber.energy_balance_func>`
        ============================  =============  =========================================================================================================================================

        Table of parameters

        ===========  =============  ==========  =======================================================================================
        Parameter    Description    Quantity    Method
        ===========  =============  ==========  =======================================================================================
        lamb         :code:`None`   ratio       :py:meth:`lambda_func <tespy.components.combustion.base.CombustionChamber.lambda_func>`
        ti           :code:`None`   heat        :py:meth:`ti_func <tespy.components.combustion.base.CombustionChamber.ti_func>`
        ===========  =============  ==========  =======================================================================================




        .. rubric:: DiabaticCombustionChamber

        Class documentation and example: :py:class:`DiabaticCombustionChamber <tespy.components.combustion.diabatic.DiabaticCombustionChamber>`

        Table of constraints

        =========================  =============  =====================================================================================================
        Parameter                  Description    Method
        =========================  =============  =====================================================================================================
        mass_flow_constraints      :code:`None`   :py:meth:`mass_flow_func <tespy.components.combustion.base.CombustionChamber.mass_flow_func>`
        stoichiometry_constraints  :code:`None`   :py:meth:`stoichiometry_func <tespy.components.combustion.base.CombustionChamber.stoichiometry_func>`
        =========================  =============  =====================================================================================================

        Table of parameters

        ===========  =============  ==========  ===================================================================================================================
        Parameter    Description    Quantity    Method
        ===========  =============  ==========  ===================================================================================================================
        lamb         :code:`None`   ratio       :py:meth:`lambda_func <tespy.components.combustion.base.CombustionChamber.lambda_func>`
        ti           :code:`None`   heat        :py:meth:`ti_func <tespy.components.combustion.base.CombustionChamber.ti_func>`
        pr           :code:`None`   ratio       :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure    :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        eta          :code:`None`   efficiency  :py:meth:`energy_balance_func <tespy.components.combustion.diabatic.DiabaticCombustionChamber.energy_balance_func>`
        Qloss        :code:`None`   heat        :code:`None`
        ===========  =============  ==========  ===================================================================================================================




        .. rubric:: CombustionEngine

        Class documentation and example: :py:class:`CombustionEngine <tespy.components.combustion.engine.CombustionEngine>`

        Table of constraints

        =============================  =============  =========================================================================================================================================
        Parameter                      Description    Method
        =============================  =============  =========================================================================================================================================
        mass_flow_constraints          :code:`None`   :py:meth:`mass_flow_func <tespy.components.combustion.base.CombustionChamber.mass_flow_func>`
        reactor_pressure_constraints   :code:`None`   :py:meth:`combustion_pressure_structure_matrix <tespy.components.combustion.base.CombustionChamber.combustion_pressure_structure_matrix>`
        stoichiometry_constraints      :code:`None`   :py:meth:`stoichiometry_func <tespy.components.combustion.base.CombustionChamber.stoichiometry_func>`
        energy_balance_constraints     :code:`None`   :py:meth:`energy_balance_func <tespy.components.combustion.engine.CombustionEngine.energy_balance_func>`
        power_constraints              :code:`None`   :py:meth:`tiP_char_func <tespy.components.combustion.engine.CombustionEngine.tiP_char_func>`
        heat1_constraints              :code:`None`   :py:meth:`Q1_char_func <tespy.components.combustion.engine.CombustionEngine.Q1_char_func>`
        heat2_constraints              :code:`None`   :py:meth:`Q2_char_func <tespy.components.combustion.engine.CombustionEngine.Q2_char_func>`
        heatloss_constraints           :code:`None`   :py:meth:`Qloss_char_func <tespy.components.combustion.engine.CombustionEngine.Qloss_char_func>`
        mass_flow_cooling_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.combustion.engine.CombustionEngine.variable_equality_structure_matrix>`
        fluid_cooling_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.combustion.engine.CombustionEngine.variable_equality_structure_matrix>`
        =============================  =============  =========================================================================================================================================

        Table of parameters

        ===========  =============  ============  =========================================================================================
        Parameter    Description    Quantity      Method
        ===========  =============  ============  =========================================================================================
        lamb         :code:`None`   ratio         :py:meth:`lambda_func <tespy.components.combustion.base.CombustionChamber.lambda_func>`
        ti           :code:`None`   heat          :py:meth:`ti_func <tespy.components.combustion.base.CombustionChamber.ti_func>`
        P            :code:`None`   power         :code:`None`
        Q1           :code:`None`   heat          :py:meth:`Q1_func <tespy.components.combustion.engine.CombustionEngine.Q1_func>`
        Q2           :code:`None`   heat          :py:meth:`Q2_func <tespy.components.combustion.engine.CombustionEngine.Q2_func>`
        Qloss        :code:`None`   heat          :code:`None`
        pr1          :code:`None`   ratio         :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        pr2          :code:`None`   ratio         :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp1          :code:`None`   pressure      :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        dp2          :code:`None`   pressure      :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta1        :code:`None`   :code:`None`  :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        zeta2        :code:`None`   :code:`None`  :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        eta_mech     :code:`None`   :code:`None`  :code:`None`
        T_v_inner    :code:`None`   :code:`None`  :code:`None`
        ===========  =============  ============  =========================================================================================

        Table of characteristic lines and maps

        ===========  =============  ============
        Parameter    Description    Method
        ===========  =============  ============
        tiP_char     :code:`None`   :code:`None`
        Q1_char      :code:`None`   :code:`None`
        Q2_char      :code:`None`   :code:`None`
        Qloss_char   :code:`None`   :code:`None`
        ===========  =============  ============



    .. tab-item:: displacementmachinery

        .. rubric:: DisplacementMachine

        Class documentation and example: :py:class:`DisplacementMachine <tespy.components.displacementmachinery.base.DisplacementMachine>`

        Table of constraints

        =====================  =============  =======================================================================================================================
        Parameter              Description    Method
        =====================  =============  =======================================================================================================================
        mass_flow_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        =====================  =============  =======================================================================================================================

        Table of parameters

        ===========  =============  ==========  ====================================================================================================================
        Parameter    Description    Quantity    Method
        ===========  =============  ==========  ====================================================================================================================
        P            :code:`None`   power       :py:meth:`energy_balance_func <tespy.components.displacementmachinery.base.DisplacementMachine.energy_balance_func>`
        pr           :code:`None`   ratio       :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure    :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        ===========  =============  ==========  ====================================================================================================================




        .. rubric:: PolynomialCompressor

        Class documentation and example: :py:class:`PolynomialCompressor <tespy.components.displacementmachinery.polynomial_compressor.PolynomialCompressor>`

        Table of constraints

        =====================  =============  =======================================================================================================================
        Parameter              Description    Method
        =====================  =============  =======================================================================================================================
        mass_flow_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        =====================  =============  =======================================================================================================================

        Table of parameters

        =================  =============  ============  =========================================================================================
        Parameter          Description    Quantity      Method
        =================  =============  ============  =========================================================================================
        P                  :code:`None`   power         :code:`None`
        pr                 :code:`None`   ratio         :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp                 :code:`None`   pressure      :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        Q_diss             :code:`None`   heat          :code:`None`
        eta_vol            :code:`None`   efficiency    :code:`None`
        dissipation_ratio  :code:`None`   ratio         :code:`None`
        Q_diss_rel         :code:`None`   ratio         :code:`None`
        rpm                :code:`None`   :code:`None`  :code:`None`
        reference_state    :code:`None`   :code:`None`  :code:`None`
        eta_s_poly         :code:`None`   :code:`None`  :code:`None`
        eta_vol_poly       :code:`None`   :code:`None`  :code:`None`
        eta_s              :code:`None`   efficiency    :code:`None`
        =================  =============  ============  =========================================================================================

        Table of parameter groups

        ====================  =============  ==========================================================  ==================================================================================================================================================
        Parameter             Description    Required parameters                                         Method
        ====================  =============  ==========================================================  ==================================================================================================================================================
        eta_vol_poly_group    :code:`None`   :code:`reference_state`, :code:`eta_vol_poly`, :code:`rpm`  :py:meth:`eta_vol_poly_group_func <tespy.components.displacementmachinery.polynomial_compressor.PolynomialCompressor.eta_vol_poly_group_func>`
        eta_vol_group         :code:`None`   :code:`reference_state`, :code:`eta_vol`, :code:`rpm`       :py:meth:`eta_vol_group_func <tespy.components.displacementmachinery.polynomial_compressor.PolynomialCompressor.eta_vol_group_func>`
        eta_s_poly_group      :code:`None`   :code:`eta_s_poly`, :code:`dissipation_ratio`               :py:meth:`eta_s_poly_group_func <tespy.components.displacementmachinery.polynomial_compressor.PolynomialCompressor.eta_s_poly_group_func>`
        eta_s_group           :code:`None`   :code:`eta_s`, :code:`dissipation_ratio`                    :py:meth:`eta_s_group_func <tespy.components.displacementmachinery.polynomial_compressor.PolynomialCompressor.eta_s_group_func>`
        energy_balance_group  :code:`None`   :code:`P`, :code:`dissipation_ratio`                        :py:meth:`energy_balance_group_func <tespy.components.displacementmachinery.polynomial_compressor.PolynomialCompressor.energy_balance_group_func>`
        ====================  =============  ==========================================================  ==================================================================================================================================================




        .. rubric:: PolynomialCompressorWithCooling

        Class documentation and example: :py:class:`PolynomialCompressorWithCooling <tespy.components.displacementmachinery.polynomial_compressor_with_cooling.PolynomialCompressorWithCooling>`

        Table of constraints

        ==================================  =============  ==============================================================================================================================================================================
        Parameter                           Description    Method
        ==================================  =============  ==============================================================================================================================================================================
        mass_flow_constraints               :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints                   :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        cooling_energy_balance_constraints  :code:`None`   :py:meth:`cooling_energy_balance_func <tespy.components.displacementmachinery.polynomial_compressor_with_cooling.PolynomialCompressorWithCooling.cooling_energy_balance_func>`
        ==================================  =============  ==============================================================================================================================================================================

        Table of parameters

        =================  =============  ======================  =========================================================================================
        Parameter          Description    Quantity                Method
        =================  =============  ======================  =========================================================================================
        P                  :code:`None`   power                   :code:`None`
        pr                 :code:`None`   ratio                   :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp                 :code:`None`   pressure                :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        Q_diss             :code:`None`   heat                    :code:`None`
        eta_vol            :code:`None`   efficiency              :code:`None`
        dissipation_ratio  :code:`None`   ratio                   :code:`None`
        Q_diss_rel         :code:`None`   ratio                   :code:`None`
        rpm                :code:`None`   :code:`None`            :code:`None`
        reference_state    :code:`None`   :code:`None`            :code:`None`
        eta_s_poly         :code:`None`   :code:`None`            :code:`None`
        eta_vol_poly       :code:`None`   :code:`None`            :code:`None`
        eta_s              :code:`None`   efficiency              :code:`None`
        eta_recovery       :code:`None`   efficiency              :code:`None`
        td_minimal         :code:`None`   temperature_difference  :code:`None`
        dp_cooling         :code:`None`   pressure                :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        pr_cooling         :code:`None`   ratio                   :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        =================  =============  ======================  =========================================================================================

        Table of parameter groups

        ====================  =============  ==========================================================  ==================================================================================================================================================
        Parameter             Description    Required parameters                                         Method
        ====================  =============  ==========================================================  ==================================================================================================================================================
        eta_vol_poly_group    :code:`None`   :code:`reference_state`, :code:`eta_vol_poly`, :code:`rpm`  :py:meth:`eta_vol_poly_group_func <tespy.components.displacementmachinery.polynomial_compressor.PolynomialCompressor.eta_vol_poly_group_func>`
        eta_vol_group         :code:`None`   :code:`reference_state`, :code:`eta_vol`, :code:`rpm`       :py:meth:`eta_vol_group_func <tespy.components.displacementmachinery.polynomial_compressor.PolynomialCompressor.eta_vol_group_func>`
        eta_s_poly_group      :code:`None`   :code:`eta_s_poly`, :code:`dissipation_ratio`               :py:meth:`eta_s_poly_group_func <tespy.components.displacementmachinery.polynomial_compressor.PolynomialCompressor.eta_s_poly_group_func>`
        eta_s_group           :code:`None`   :code:`eta_s`, :code:`dissipation_ratio`                    :py:meth:`eta_s_group_func <tespy.components.displacementmachinery.polynomial_compressor.PolynomialCompressor.eta_s_group_func>`
        energy_balance_group  :code:`None`   :code:`P`, :code:`dissipation_ratio`                        :py:meth:`energy_balance_group_func <tespy.components.displacementmachinery.polynomial_compressor.PolynomialCompressor.energy_balance_group_func>`
        ====================  =============  ==========================================================  ==================================================================================================================================================



    .. tab-item:: heat_exchangers

        .. rubric:: HeatExchanger

        Class documentation and example: :py:class:`HeatExchanger <tespy.components.heat_exchangers.base.HeatExchanger>`

        Table of constraints

        ==========================  =============  =======================================================================================================================
        Parameter                   Description    Method
        ==========================  =============  =======================================================================================================================
        mass_flow_constraints       :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints           :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        energy_balance_constraints  :code:`None`   :py:meth:`energy_balance_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_func>`
        ==========================  =============  =======================================================================================================================

        Table of parameters

        ===========  =============  =========================  ================================================================================================================
        Parameter    Description    Quantity                   Method
        ===========  =============  =========================  ================================================================================================================
        Q            :code:`None`   heat                       :py:meth:`energy_balance_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func>`
        kA           :code:`None`   heat_transfer_coefficient  :py:meth:`kA_func <tespy.components.heat_exchangers.base.HeatExchanger.kA_func>`
        td_log       :code:`None`   temperature_difference     :code:`None`
        ttd_u        :code:`None`   temperature_difference     :py:meth:`ttd_u_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_u_func>`
        ttd_l        :code:`None`   temperature_difference     :py:meth:`ttd_l_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_l_func>`
        ttd_min      :code:`None`   temperature_difference     :py:meth:`ttd_min_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_min_func>`
        pr1          :code:`None`   ratio                      :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        pr2          :code:`None`   ratio                      :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp1          :code:`None`   pressure                   :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        dp2          :code:`None`   pressure                   :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta1        :code:`None`   :code:`None`               :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        zeta2        :code:`None`   :code:`None`               :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        eff_cold     :code:`None`   efficiency                 :py:meth:`eff_cold_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_cold_func>`
        eff_hot      :code:`None`   efficiency                 :py:meth:`eff_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_hot_func>`
        eff_max      :code:`None`   efficiency                 :py:meth:`eff_max_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_max_func>`
        ===========  =============  =========================  ================================================================================================================

        Table of parameter groups

        ===========  =============  ==================================  ==========================================================================================
        Parameter    Description    Required parameters                 Method
        ===========  =============  ==================================  ==========================================================================================
        kA_char      :code:`None`   :code:`kA_char1`, :code:`kA_char2`  :py:meth:`kA_char_func <tespy.components.heat_exchangers.base.HeatExchanger.kA_char_func>`
        ===========  =============  ==================================  ==========================================================================================

        Table of characteristic lines and maps

        ===========  =============  ============
        Parameter    Description    Method
        ===========  =============  ============
        kA_char1     :code:`None`   :code:`None`
        kA_char2     :code:`None`   :code:`None`
        ===========  =============  ============




        .. rubric:: Condenser

        Class documentation and example: :py:class:`Condenser <tespy.components.heat_exchangers.condenser.Condenser>`

        Table of constraints

        ==========================  =============  =======================================================================================================================
        Parameter                   Description    Method
        ==========================  =============  =======================================================================================================================
        mass_flow_constraints       :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints           :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        energy_balance_constraints  :code:`None`   :py:meth:`energy_balance_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_func>`
        ==========================  =============  =======================================================================================================================

        Table of parameters

        ===========  =============  =========================  ================================================================================================================
        Parameter    Description    Quantity                   Method
        ===========  =============  =========================  ================================================================================================================
        Q            :code:`None`   heat                       :py:meth:`energy_balance_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func>`
        kA           :code:`None`   heat_transfer_coefficient  :py:meth:`kA_func <tespy.components.heat_exchangers.base.HeatExchanger.kA_func>`
        td_log       :code:`None`   temperature_difference     :code:`None`
        ttd_u        :code:`None`   temperature_difference     :py:meth:`ttd_u_func <tespy.components.heat_exchangers.condenser.Condenser.ttd_u_func>`
        ttd_l        :code:`None`   temperature_difference     :py:meth:`ttd_l_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_l_func>`
        ttd_min      :code:`None`   temperature_difference     :py:meth:`ttd_min_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_min_func>`
        pr1          :code:`None`   ratio                      :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        pr2          :code:`None`   ratio                      :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp1          :code:`None`   pressure                   :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        dp2          :code:`None`   pressure                   :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta1        :code:`None`   :code:`None`               :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        zeta2        :code:`None`   :code:`None`               :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        eff_cold     :code:`None`   efficiency                 :py:meth:`eff_cold_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_cold_func>`
        eff_hot      :code:`None`   efficiency                 :py:meth:`eff_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_hot_func>`
        eff_max      :code:`None`   efficiency                 :py:meth:`eff_max_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_max_func>`
        subcooling   :code:`None`   :code:`None`               :py:meth:`subcooling_func <tespy.components.heat_exchangers.condenser.Condenser.subcooling_func>`
        ===========  =============  =========================  ================================================================================================================

        Table of parameter groups

        ===========  =============  ==================================  ===========================================================================================
        Parameter    Description    Required parameters                 Method
        ===========  =============  ==================================  ===========================================================================================
        kA_char      :code:`None`   :code:`kA_char1`, :code:`kA_char2`  :py:meth:`kA_char_func <tespy.components.heat_exchangers.condenser.Condenser.kA_char_func>`
        ===========  =============  ==================================  ===========================================================================================

        Table of characteristic lines and maps

        ===========  =============  ============
        Parameter    Description    Method
        ===========  =============  ============
        kA_char1     :code:`None`   :code:`None`
        kA_char2     :code:`None`   :code:`None`
        ===========  =============  ============




        .. rubric:: Desuperheater

        Class documentation and example: :py:class:`Desuperheater <tespy.components.heat_exchangers.desuperheater.Desuperheater>`

        Table of constraints

        ==========================  =============  =======================================================================================================================
        Parameter                   Description    Method
        ==========================  =============  =======================================================================================================================
        mass_flow_constraints       :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints           :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        energy_balance_constraints  :code:`None`   :py:meth:`energy_balance_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_func>`
        saturated_gas_constraints   :code:`None`   :py:meth:`saturated_gas_func <tespy.components.heat_exchangers.desuperheater.Desuperheater.saturated_gas_func>`
        ==========================  =============  =======================================================================================================================

        Table of parameters

        ===========  =============  =========================  ================================================================================================================
        Parameter    Description    Quantity                   Method
        ===========  =============  =========================  ================================================================================================================
        Q            :code:`None`   heat                       :py:meth:`energy_balance_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func>`
        kA           :code:`None`   heat_transfer_coefficient  :py:meth:`kA_func <tespy.components.heat_exchangers.base.HeatExchanger.kA_func>`
        td_log       :code:`None`   temperature_difference     :code:`None`
        ttd_u        :code:`None`   temperature_difference     :py:meth:`ttd_u_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_u_func>`
        ttd_l        :code:`None`   temperature_difference     :py:meth:`ttd_l_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_l_func>`
        ttd_min      :code:`None`   temperature_difference     :py:meth:`ttd_min_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_min_func>`
        pr1          :code:`None`   ratio                      :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        pr2          :code:`None`   ratio                      :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp1          :code:`None`   pressure                   :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        dp2          :code:`None`   pressure                   :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta1        :code:`None`   :code:`None`               :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        zeta2        :code:`None`   :code:`None`               :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        eff_cold     :code:`None`   efficiency                 :py:meth:`eff_cold_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_cold_func>`
        eff_hot      :code:`None`   efficiency                 :py:meth:`eff_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_hot_func>`
        eff_max      :code:`None`   efficiency                 :py:meth:`eff_max_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_max_func>`
        ===========  =============  =========================  ================================================================================================================

        Table of parameter groups

        ===========  =============  ==================================  ==========================================================================================
        Parameter    Description    Required parameters                 Method
        ===========  =============  ==================================  ==========================================================================================
        kA_char      :code:`None`   :code:`kA_char1`, :code:`kA_char2`  :py:meth:`kA_char_func <tespy.components.heat_exchangers.base.HeatExchanger.kA_char_func>`
        ===========  =============  ==================================  ==========================================================================================

        Table of characteristic lines and maps

        ===========  =============  ============
        Parameter    Description    Method
        ===========  =============  ============
        kA_char1     :code:`None`   :code:`None`
        kA_char2     :code:`None`   :code:`None`
        ===========  =============  ============




        .. rubric:: SectionedHeatExchanger

        Class documentation and example: :py:class:`SectionedHeatExchanger <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger>`

        Table of constraints

        ==========================  =============  =======================================================================================================================
        Parameter                   Description    Method
        ==========================  =============  =======================================================================================================================
        mass_flow_constraints       :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints           :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        energy_balance_constraints  :code:`None`   :py:meth:`energy_balance_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_func>`
        ==========================  =============  =======================================================================================================================

        Table of parameters

        =================  =============  =========================  ================================================================================================================
        Parameter          Description    Quantity                   Method
        =================  =============  =========================  ================================================================================================================
        Q                  :code:`None`   heat                       :py:meth:`energy_balance_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func>`
        kA                 :code:`None`   heat_transfer_coefficient  :py:meth:`kA_func <tespy.components.heat_exchangers.base.HeatExchanger.kA_func>`
        td_log             :code:`None`   temperature_difference     :code:`None`
        ttd_u              :code:`None`   temperature_difference     :py:meth:`ttd_u_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_u_func>`
        ttd_l              :code:`None`   temperature_difference     :py:meth:`ttd_l_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_l_func>`
        ttd_min            :code:`None`   temperature_difference     :py:meth:`ttd_min_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_min_func>`
        pr1                :code:`None`   ratio                      :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        pr2                :code:`None`   ratio                      :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp1                :code:`None`   pressure                   :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        dp2                :code:`None`   pressure                   :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta1              :code:`None`   :code:`None`               :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        zeta2              :code:`None`   :code:`None`               :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        eff_cold           :code:`None`   efficiency                 :py:meth:`eff_cold_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_cold_func>`
        eff_hot            :code:`None`   efficiency                 :py:meth:`eff_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_hot_func>`
        eff_max            :code:`None`   efficiency                 :py:meth:`eff_max_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_max_func>`
        num_sections       :code:`None`   :code:`None`               :code:`None`
        UA                 :code:`None`   heat_transfer_coefficient  :py:meth:`UA_func <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.UA_func>`
        refrigerant_index  :code:`None`   :code:`None`               :code:`None`
        re_exp_r           :code:`None`   :code:`None`               :code:`None`
        re_exp_sf          :code:`None`   :code:`None`               :code:`None`
        alpha_ratio        :code:`None`   ratio                      :code:`None`
        area_ratio         :code:`None`   ratio                      :code:`None`
        td_pinch           :code:`None`   temperature_difference     :py:meth:`td_pinch_func <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.td_pinch_func>`
        =================  =============  =========================  ================================================================================================================

        Table of parameter groups

        =============  =============  ============================================================================  ====================================================================================================================
        Parameter      Description    Required parameters                                                           Method
        =============  =============  ============================================================================  ====================================================================================================================
        kA_char        :code:`None`   :code:`kA_char1`, :code:`kA_char2`                                            :py:meth:`kA_char_func <tespy.components.heat_exchangers.base.HeatExchanger.kA_char_func>`
        UA_cecchinato  :code:`None`   :code:`re_exp_r`, :code:`re_exp_sf`, :code:`alpha_ratio`, :code:`area_ratio`  :py:meth:`UA_cecchinato_func <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.UA_cecchinato_func>`
        =============  =============  ============================================================================  ====================================================================================================================

        Table of characteristic lines and maps

        ===========  =============  ============
        Parameter    Description    Method
        ===========  =============  ============
        kA_char1     :code:`None`   :code:`None`
        kA_char2     :code:`None`   :code:`None`
        ===========  =============  ============




        .. rubric:: MovingBoundaryHeatExchanger

        Class documentation and example: :py:class:`MovingBoundaryHeatExchanger <tespy.components.heat_exchangers.movingboundary.MovingBoundaryHeatExchanger>`

        Table of constraints

        ==========================  =============  =======================================================================================================================
        Parameter                   Description    Method
        ==========================  =============  =======================================================================================================================
        mass_flow_constraints       :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints           :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        energy_balance_constraints  :code:`None`   :py:meth:`energy_balance_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_func>`
        ==========================  =============  =======================================================================================================================

        Table of parameters

        =================  =============  =========================  ================================================================================================================
        Parameter          Description    Quantity                   Method
        =================  =============  =========================  ================================================================================================================
        Q                  :code:`None`   heat                       :py:meth:`energy_balance_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func>`
        kA                 :code:`None`   heat_transfer_coefficient  :py:meth:`kA_func <tespy.components.heat_exchangers.base.HeatExchanger.kA_func>`
        td_log             :code:`None`   temperature_difference     :code:`None`
        ttd_u              :code:`None`   temperature_difference     :py:meth:`ttd_u_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_u_func>`
        ttd_l              :code:`None`   temperature_difference     :py:meth:`ttd_l_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_l_func>`
        ttd_min            :code:`None`   temperature_difference     :py:meth:`ttd_min_func <tespy.components.heat_exchangers.base.HeatExchanger.ttd_min_func>`
        pr1                :code:`None`   ratio                      :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        pr2                :code:`None`   ratio                      :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp1                :code:`None`   pressure                   :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        dp2                :code:`None`   pressure                   :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta1              :code:`None`   :code:`None`               :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        zeta2              :code:`None`   :code:`None`               :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        eff_cold           :code:`None`   efficiency                 :py:meth:`eff_cold_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_cold_func>`
        eff_hot            :code:`None`   efficiency                 :py:meth:`eff_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_hot_func>`
        eff_max            :code:`None`   efficiency                 :py:meth:`eff_max_func <tespy.components.heat_exchangers.base.HeatExchanger.eff_max_func>`
        UA                 :code:`None`   heat_transfer_coefficient  :py:meth:`UA_func <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.UA_func>`
        refrigerant_index  :code:`None`   :code:`None`               :code:`None`
        re_exp_r           :code:`None`   :code:`None`               :code:`None`
        re_exp_sf          :code:`None`   :code:`None`               :code:`None`
        alpha_ratio        :code:`None`   ratio                      :code:`None`
        area_ratio         :code:`None`   ratio                      :code:`None`
        td_pinch           :code:`None`   temperature_difference     :py:meth:`td_pinch_func <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.td_pinch_func>`
        =================  =============  =========================  ================================================================================================================

        Table of parameter groups

        =============  =============  ============================================================================  ====================================================================================================================
        Parameter      Description    Required parameters                                                           Method
        =============  =============  ============================================================================  ====================================================================================================================
        kA_char        :code:`None`   :code:`kA_char1`, :code:`kA_char2`                                            :py:meth:`kA_char_func <tespy.components.heat_exchangers.base.HeatExchanger.kA_char_func>`
        UA_cecchinato  :code:`None`   :code:`re_exp_r`, :code:`re_exp_sf`, :code:`alpha_ratio`, :code:`area_ratio`  :py:meth:`UA_cecchinato_func <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.UA_cecchinato_func>`
        =============  =============  ============================================================================  ====================================================================================================================

        Table of characteristic lines and maps

        ===========  =============  ============
        Parameter    Description    Method
        ===========  =============  ============
        kA_char1     :code:`None`   :code:`None`
        kA_char2     :code:`None`   :code:`None`
        ===========  =============  ============




        .. rubric:: SimpleHeatExchanger

        Class documentation and example: :py:class:`SimpleHeatExchanger <tespy.components.heat_exchangers.simple.SimpleHeatExchanger>`

        Table of constraints

        =====================  =============  =======================================================================================================================
        Parameter              Description    Method
        =====================  =============  =======================================================================================================================
        mass_flow_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        =====================  =============  =======================================================================================================================

        Table of parameters

        ========================  =============  =========================  ================================================================================================================
        Parameter                 Description    Quantity                   Method
        ========================  =============  =========================  ================================================================================================================
        power_connector_location  :code:`None`   :code:`None`               :code:`None`
        Q                         :code:`None`   heat                       :py:meth:`energy_balance_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_balance_func>`
        pr                        :code:`None`   ratio                      :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp                        :code:`None`   pressure                   :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta                      :code:`None`   :code:`None`               :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        D                         :code:`None`   length                     :code:`None`
        L                         :code:`None`   length                     :code:`None`
        ks                        :code:`None`   length                     :code:`None`
        ks_HW                     :code:`None`   :code:`None`               :code:`None`
        kA                        :code:`None`   heat_transfer_coefficient  :code:`None`
        Tamb                      :code:`None`   temperature                :code:`None`
        dissipative               :code:`None`   :code:`None`               :code:`None`
        ========================  =============  =========================  ================================================================================================================

        Table of parameter groups

        =============  =============  ===================================  ================================================================================================================
        Parameter      Description    Required parameters                  Method
        =============  =============  ===================================  ================================================================================================================
        darcy_group    :code:`None`   :code:`L`, :code:`ks`, :code:`D`     :py:meth:`darcy_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.darcy_func>`
        hw_group       :code:`None`   :code:`L`, :code:`ks_HW`, :code:`D`  :py:meth:`hazen_williams_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.hazen_williams_func>`
        kA_group       :code:`None`   :code:`kA`, :code:`Tamb`             :py:meth:`kA_group_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.kA_group_func>`
        kA_char_group  :code:`None`   :code:`kA_char`, :code:`Tamb`        :py:meth:`kA_char_group_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.kA_char_group_func>`
        =============  =============  ===================================  ================================================================================================================

        Table of characteristic lines and maps

        ===========  =============  ============
        Parameter    Description    Method
        ===========  =============  ============
        kA_char      :code:`None`   :code:`None`
        ===========  =============  ============




        .. rubric:: ParabolicTrough

        Class documentation and example: :py:class:`ParabolicTrough <tespy.components.heat_exchangers.parabolic_trough.ParabolicTrough>`

        Table of constraints

        =====================  =============  =======================================================================================================================
        Parameter              Description    Method
        =====================  =============  =======================================================================================================================
        mass_flow_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        =====================  =============  =======================================================================================================================

        Table of parameters

        ========================  =============  ============  ================================================================================================================
        Parameter                 Description    Quantity      Method
        ========================  =============  ============  ================================================================================================================
        power_connector_location  :code:`None`   :code:`None`  :code:`None`
        Q                         :code:`None`   heat          :py:meth:`energy_balance_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_balance_func>`
        pr                        :code:`None`   ratio         :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp                        :code:`None`   pressure      :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta                      :code:`None`   :code:`None`  :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        D                         :code:`None`   length        :code:`None`
        L                         :code:`None`   length        :code:`None`
        ks                        :code:`None`   length        :code:`None`
        ks_HW                     :code:`None`   :code:`None`  :code:`None`
        Tamb                      :code:`None`   temperature   :code:`None`
        dissipative               :code:`None`   :code:`None`  :code:`None`
        E                         :code:`None`   heat          :code:`None`
        A                         :code:`None`   area          :code:`None`
        eta_opt                   :code:`None`   efficiency    :code:`None`
        c_1                       :code:`None`   :code:`None`  :code:`None`
        c_2                       :code:`None`   :code:`None`  :code:`None`
        iam_1                     :code:`None`   :code:`None`  :code:`None`
        iam_2                     :code:`None`   :code:`None`  :code:`None`
        aoi                       :code:`None`   angle         :code:`None`
        doc                       :code:`None`   ratio         :code:`None`
        Q_loss                    :code:`None`   heat          :code:`None`
        ========================  =============  ============  ================================================================================================================

        Table of parameter groups

        ============  =============  =====================================================================================================================================  ==================================================================================================================
        Parameter     Description    Required parameters                                                                                                                    Method
        ============  =============  =====================================================================================================================================  ==================================================================================================================
        darcy_group   :code:`None`   :code:`L`, :code:`ks`, :code:`D`                                                                                                       :py:meth:`darcy_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.darcy_func>`
        hw_group      :code:`None`   :code:`L`, :code:`ks_HW`, :code:`D`                                                                                                    :py:meth:`hazen_williams_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.hazen_williams_func>`
        energy_group  :code:`None`   :code:`E`, :code:`eta_opt`, :code:`aoi`, :code:`doc`, :code:`c_1`, :code:`c_2`, :code:`iam_1`, :code:`iam_2`, :code:`A`, :code:`Tamb`  :py:meth:`energy_group_func <tespy.components.heat_exchangers.parabolic_trough.ParabolicTrough.energy_group_func>`
        ============  =============  =====================================================================================================================================  ==================================================================================================================




        .. rubric:: ParallelFlowHeatExchanger

        Class documentation and example: :py:class:`ParallelFlowHeatExchanger <tespy.components.heat_exchangers.parallel.ParallelFlowHeatExchanger>`

        Table of constraints

        ==========================  =============  =======================================================================================================================
        Parameter                   Description    Method
        ==========================  =============  =======================================================================================================================
        mass_flow_constraints       :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints           :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        energy_balance_constraints  :code:`None`   :py:meth:`energy_balance_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_func>`
        ==========================  =============  =======================================================================================================================

        Table of parameters

        ===========  =============  =========================  ================================================================================================================
        Parameter    Description    Quantity                   Method
        ===========  =============  =========================  ================================================================================================================
        Q            :code:`None`   heat                       :py:meth:`energy_balance_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func>`
        kA           :code:`None`   heat_transfer_coefficient  :py:meth:`kA_func <tespy.components.heat_exchangers.base.HeatExchanger.kA_func>`
        td_log       :code:`None`   temperature_difference     :code:`None`
        ttd_u        :code:`None`   temperature_difference     :py:meth:`ttd_u_func <tespy.components.heat_exchangers.parallel.ParallelFlowHeatExchanger.ttd_u_func>`
        ttd_l        :code:`None`   temperature_difference     :py:meth:`ttd_l_func <tespy.components.heat_exchangers.parallel.ParallelFlowHeatExchanger.ttd_l_func>`
        pr1          :code:`None`   ratio                      :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        pr2          :code:`None`   ratio                      :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp1          :code:`None`   pressure                   :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        dp2          :code:`None`   pressure                   :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta1        :code:`None`   :code:`None`               :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        zeta2        :code:`None`   :code:`None`               :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        ===========  =============  =========================  ================================================================================================================

        Table of parameter groups

        ===========  =============  ==================================  ==========================================================================================
        Parameter    Description    Required parameters                 Method
        ===========  =============  ==================================  ==========================================================================================
        kA_char      :code:`None`   :code:`kA_char1`, :code:`kA_char2`  :py:meth:`kA_char_func <tespy.components.heat_exchangers.base.HeatExchanger.kA_char_func>`
        ===========  =============  ==================================  ==========================================================================================

        Table of characteristic lines and maps

        ===========  =============  ============
        Parameter    Description    Method
        ===========  =============  ============
        kA_char1     :code:`None`   :code:`None`
        kA_char2     :code:`None`   :code:`None`
        ===========  =============  ============




        .. rubric:: SolarCollector

        Class documentation and example: :py:class:`SolarCollector <tespy.components.heat_exchangers.solar_collector.SolarCollector>`

        Table of constraints

        =====================  =============  =======================================================================================================================
        Parameter              Description    Method
        =====================  =============  =======================================================================================================================
        mass_flow_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        =====================  =============  =======================================================================================================================

        Table of parameters

        ========================  =============  ============  ================================================================================================================
        Parameter                 Description    Quantity      Method
        ========================  =============  ============  ================================================================================================================
        power_connector_location  :code:`None`   :code:`None`  :code:`None`
        Q                         :code:`None`   heat          :py:meth:`energy_balance_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_balance_func>`
        pr                        :code:`None`   ratio         :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp                        :code:`None`   pressure      :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta                      :code:`None`   :code:`None`  :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        D                         :code:`None`   length        :code:`None`
        L                         :code:`None`   length        :code:`None`
        ks                        :code:`None`   length        :code:`None`
        ks_HW                     :code:`None`   :code:`None`  :code:`None`
        Tamb                      :code:`None`   temperature   :code:`None`
        dissipative               :code:`None`   :code:`None`  :code:`None`
        E                         :code:`None`   heat          :code:`None`
        A                         :code:`None`   area          :code:`None`
        eta_opt                   :code:`None`   efficiency    :code:`None`
        lkf_lin                   :code:`None`   :code:`None`  :code:`None`
        lkf_quad                  :code:`None`   :code:`None`  :code:`None`
        Q_loss                    :code:`None`   temperature   :code:`None`
        ========================  =============  ============  ================================================================================================================

        Table of parameter groups

        ============  =============  ======================================================================================  ================================================================================================================
        Parameter     Description    Required parameters                                                                     Method
        ============  =============  ======================================================================================  ================================================================================================================
        darcy_group   :code:`None`   :code:`L`, :code:`ks`, :code:`D`                                                        :py:meth:`darcy_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.darcy_func>`
        hw_group      :code:`None`   :code:`L`, :code:`ks_HW`, :code:`D`                                                     :py:meth:`hazen_williams_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.hazen_williams_func>`
        energy_group  :code:`None`   :code:`E`, :code:`eta_opt`, :code:`lkf_lin`, :code:`lkf_quad`, :code:`A`, :code:`Tamb`  :py:meth:`energy_group_func <tespy.components.heat_exchangers.solar_collector.SolarCollector.energy_group_func>`
        ============  =============  ======================================================================================  ================================================================================================================



    .. tab-item:: nodes

        .. rubric:: NodeBase

        Class documentation and example: :py:class:`NodeBase <tespy.components.nodes.base.NodeBase>`

        Table of constraints

        =====================  =============  =======================================================================================================================
        Parameter              Description    Method
        =====================  =============  =======================================================================================================================
        mass_flow_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        =====================  =============  =======================================================================================================================




        .. rubric:: DropletSeparator

        Class documentation and example: :py:class:`DropletSeparator <tespy.components.nodes.droplet_separator.DropletSeparator>`

        Table of constraints

        ==========================  =============  ====================================================================================================================
        Parameter                   Description    Method
        ==========================  =============  ====================================================================================================================
        mass_flow_constraints       :code:`None`   :py:meth:`mass_flow_func <tespy.components.nodes.base.NodeBase.mass_flow_func>`
        energy_balance_constraints  :code:`None`   :py:meth:`energy_balance_func <tespy.components.nodes.droplet_separator.DropletSeparator.energy_balance_func>`
        pressure_constraints        :code:`None`   :py:meth:`pressure_structure_matrix <tespy.components.nodes.base.NodeBase.pressure_structure_matrix>`
        outlet_constraint_liquid    :code:`None`   :py:meth:`saturated_outlet_func <tespy.components.nodes.droplet_separator.DropletSeparator.saturated_outlet_func>`
        outlet_constraint_gas       :code:`None`   :py:meth:`saturated_outlet_func <tespy.components.nodes.droplet_separator.DropletSeparator.saturated_outlet_func>`
        fluid_constraints           :code:`None`   :py:meth:`fluid_structure_matrix <tespy.components.nodes.droplet_separator.DropletSeparator.fluid_structure_matrix>`
        ==========================  =============  ====================================================================================================================




        .. rubric:: Drum

        Class documentation and example: :py:class:`Drum <tespy.components.nodes.drum.Drum>`

        Table of constraints

        ==========================  =============  ====================================================================================================================
        Parameter                   Description    Method
        ==========================  =============  ====================================================================================================================
        mass_flow_constraints       :code:`None`   :py:meth:`mass_flow_func <tespy.components.nodes.base.NodeBase.mass_flow_func>`
        energy_balance_constraints  :code:`None`   :py:meth:`energy_balance_func <tespy.components.nodes.droplet_separator.DropletSeparator.energy_balance_func>`
        pressure_constraints        :code:`None`   :py:meth:`pressure_structure_matrix <tespy.components.nodes.base.NodeBase.pressure_structure_matrix>`
        outlet_constraint_liquid    :code:`None`   :py:meth:`saturated_outlet_func <tespy.components.nodes.droplet_separator.DropletSeparator.saturated_outlet_func>`
        outlet_constraint_gas       :code:`None`   :py:meth:`saturated_outlet_func <tespy.components.nodes.droplet_separator.DropletSeparator.saturated_outlet_func>`
        fluid_constraints           :code:`None`   :py:meth:`fluid_structure_matrix <tespy.components.nodes.droplet_separator.DropletSeparator.fluid_structure_matrix>`
        ==========================  =============  ====================================================================================================================




        .. rubric:: Merge

        Class documentation and example: :py:class:`Merge <tespy.components.nodes.merge.Merge>`

        Table of constraints

        ==========================  =============  =====================================================================================================
        Parameter                   Description    Method
        ==========================  =============  =====================================================================================================
        mass_flow_constraints       :code:`None`   :py:meth:`mass_flow_func <tespy.components.nodes.base.NodeBase.mass_flow_func>`
        fluid_constraints           :code:`None`   :py:meth:`fluid_func <tespy.components.nodes.merge.Merge.fluid_func>`
        energy_balance_constraints  :code:`None`   :py:meth:`energy_balance_func <tespy.components.nodes.merge.Merge.energy_balance_func>`
        pressure_constraints        :code:`None`   :py:meth:`pressure_structure_matrix <tespy.components.nodes.base.NodeBase.pressure_structure_matrix>`
        ==========================  =============  =====================================================================================================

        Table of parameters

        ===========  =============  ============  ============
        Parameter    Description    Quantity      Method
        ===========  =============  ============  ============
        num_in       :code:`None`   :code:`None`  :code:`None`
        ===========  =============  ============  ============




        .. rubric:: Splitter

        Class documentation and example: :py:class:`Splitter <tespy.components.nodes.splitter.Splitter>`

        Table of constraints

        ==========================  =============  =========================================================================================================
        Parameter                   Description    Method
        ==========================  =============  =========================================================================================================
        mass_flow_constraints       :code:`None`   :py:meth:`mass_flow_func <tespy.components.nodes.base.NodeBase.mass_flow_func>`
        energy_balance_constraints  :code:`None`   :py:meth:`enthalpy_structure_matrix <tespy.components.nodes.splitter.Splitter.enthalpy_structure_matrix>`
        pressure_constraints        :code:`None`   :py:meth:`pressure_structure_matrix <tespy.components.nodes.base.NodeBase.pressure_structure_matrix>`
        fluid_constraints           :code:`None`   :py:meth:`fluid_structure_matrix <tespy.components.nodes.splitter.Splitter.fluid_structure_matrix>`
        ==========================  =============  =========================================================================================================

        Table of parameters

        ===========  =============  ============  ============
        Parameter    Description    Quantity      Method
        ===========  =============  ============  ============
        num_out      :code:`None`   :code:`None`  :code:`None`
        ===========  =============  ============  ============




        .. rubric:: Node

        Class documentation and example: :py:class:`Node <tespy.components.nodes.node.Node>`

        Table of constraints

        ===========================  =============  =====================================================================================================
        Parameter                    Description    Method
        ===========================  =============  =====================================================================================================
        mass_flow_constraints        :code:`None`   :py:meth:`mass_flow_func <tespy.components.nodes.base.NodeBase.mass_flow_func>`
        pressure_constraints         :code:`None`   :py:meth:`pressure_structure_matrix <tespy.components.nodes.base.NodeBase.pressure_structure_matrix>`
        outlet_enthalpy_constraints  :code:`None`   :py:meth:`enthalpy_structure_matrix <tespy.components.nodes.node.Node.enthalpy_structure_matrix>`
        outlet_fluid_constraints     :code:`None`   :py:meth:`fluid_structure_matrix <tespy.components.nodes.node.Node.fluid_structure_matrix>`
        fluid_constraints            :code:`None`   :py:meth:`fluid_func <tespy.components.nodes.merge.Merge.fluid_func>`
        energy_balance_constraints   :code:`None`   :py:meth:`energy_balance_func <tespy.components.nodes.merge.Merge.energy_balance_func>`
        ===========================  =============  =====================================================================================================

        Table of parameters

        ===========  =============  ============  ============
        Parameter    Description    Quantity      Method
        ===========  =============  ============  ============
        num_out      :code:`None`   :code:`None`  :code:`None`
        num_in       :code:`None`   :code:`None`  :code:`None`
        ===========  =============  ============  ============




        .. rubric:: Separator

        Class documentation and example: :py:class:`Separator <tespy.components.nodes.separator.Separator>`

        Table of constraints

        ==========================  =============  =====================================================================================================
        Parameter                   Description    Method
        ==========================  =============  =====================================================================================================
        mass_flow_constraints       :code:`None`   :py:meth:`mass_flow_func <tespy.components.nodes.base.NodeBase.mass_flow_func>`
        fluid_constraints           :code:`None`   :py:meth:`fluid_func <tespy.components.nodes.separator.Separator.fluid_func>`
        energy_balance_constraints  :code:`None`   :py:meth:`energy_balance_func <tespy.components.nodes.separator.Separator.energy_balance_func>`
        pressure_constraints        :code:`None`   :py:meth:`pressure_structure_matrix <tespy.components.nodes.base.NodeBase.pressure_structure_matrix>`
        ==========================  =============  =====================================================================================================

        Table of parameters

        ===========  =============  ============  ============
        Parameter    Description    Quantity      Method
        ===========  =============  ============  ============
        num_out      :code:`None`   :code:`None`  :code:`None`
        ===========  =============  ============  ============



    .. tab-item:: piping

        .. rubric:: Pipe

        Class documentation and example: :py:class:`Pipe <tespy.components.piping.pipe.Pipe>`

        Table of constraints

        =====================  =============  =======================================================================================================================
        Parameter              Description    Method
        =====================  =============  =======================================================================================================================
        mass_flow_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        =====================  =============  =======================================================================================================================

        Table of parameters

        ========================  =============  =========================  ================================================================================================================
        Parameter                 Description    Quantity                   Method
        ========================  =============  =========================  ================================================================================================================
        power_connector_location  :code:`None`   :code:`None`               :code:`None`
        Q                         :code:`None`   heat                       :py:meth:`energy_balance_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_balance_func>`
        pr                        :code:`None`   ratio                      :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp                        :code:`None`   pressure                   :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta                      :code:`None`   :code:`None`               :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        D                         :code:`None`   length                     :code:`None`
        L                         :code:`None`   length                     :code:`None`
        ks                        :code:`None`   length                     :code:`None`
        ks_HW                     :code:`None`   :code:`None`               :code:`None`
        kA                        :code:`None`   heat_transfer_coefficient  :code:`None`
        Tamb                      :code:`None`   temperature                :code:`None`
        dissipative               :code:`None`   :code:`None`               :code:`None`
        insulation_thickness      :code:`None`   length                     :code:`None`
        insulation_tc             :code:`None`   thermal_conductivity       :code:`None`
        material                  :code:`None`   :code:`None`               :code:`None`
        pipe_thickness            :code:`None`   length                     :code:`None`
        environment_media         :code:`None`   :code:`None`               :code:`None`
        wind_velocity             :code:`None`   speed                      :code:`None`
        pipe_depth                :code:`None`   length                     :code:`None`
        ========================  =============  =========================  ================================================================================================================

        Table of parameter groups

        ======================  =============  =============================================================================================================================================================  ================================================================================================================
        Parameter               Description    Required parameters                                                                                                                                            Method
        ======================  =============  =============================================================================================================================================================  ================================================================================================================
        darcy_group             :code:`None`   :code:`L`, :code:`ks`, :code:`D`                                                                                                                               :py:meth:`darcy_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.darcy_func>`
        hw_group                :code:`None`   :code:`L`, :code:`ks_HW`, :code:`D`                                                                                                                            :py:meth:`hazen_williams_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.hazen_williams_func>`
        kA_group                :code:`None`   :code:`kA`, :code:`Tamb`                                                                                                                                       :py:meth:`kA_group_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.kA_group_func>`
        kA_char_group           :code:`None`   :code:`kA_char`, :code:`Tamb`                                                                                                                                  :py:meth:`kA_char_group_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.kA_char_group_func>`
        Q_ohc_group_surface     :code:`None`   :code:`insulation_thickness`, :code:`insulation_tc`, :code:`Tamb`, :code:`material`, :code:`pipe_thickness`, :code:`environment_media`, :code:`wind_velocity`  :py:meth:`ohc_surface_group_func <tespy.components.piping.pipe.Pipe.ohc_surface_group_func>`
        Q_ohc_group_subsurface  :code:`None`   :code:`insulation_thickness`, :code:`insulation_tc`, :code:`Tamb`, :code:`material`, :code:`pipe_thickness`, :code:`environment_media`, :code:`pipe_depth`     :py:meth:`ohc_subsurface_group_func <tespy.components.piping.pipe.Pipe.ohc_subsurface_group_func>`
        ======================  =============  =============================================================================================================================================================  ================================================================================================================

        Table of characteristic lines and maps

        ===========  =============  ============
        Parameter    Description    Method
        ===========  =============  ============
        kA_char      :code:`None`   :code:`None`
        ===========  =============  ============




        .. rubric:: Valve

        Class documentation and example: :py:class:`Valve <tespy.components.piping.valve.Valve>`

        Table of constraints

        =====================  =============  =======================================================================================================================
        Parameter              Description    Method
        =====================  =============  =======================================================================================================================
        mass_flow_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        enthalpy_constraints   :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        =====================  =============  =======================================================================================================================

        Table of parameters

        ===========  =============  ============  =========================================================================================
        Parameter    Description    Quantity      Method
        ===========  =============  ============  =========================================================================================
        pr           :code:`None`   ratio         :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure      :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta         :code:`None`   :code:`None`  :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        ===========  =============  ============  =========================================================================================

        Table of characteristic lines and maps

        ===========  =============  ==========================================================================
        Parameter    Description    Method
        ===========  =============  ==========================================================================
        dp_char      :code:`None`   :py:meth:`dp_char_func <tespy.components.piping.valve.Valve.dp_char_func>`
        ===========  =============  ==========================================================================



    .. tab-item:: power

        .. rubric:: PowerBus

        Class documentation and example: :py:class:`PowerBus <tespy.components.power.bus.PowerBus>`

        Table of constraints

        =========================  =============  ========================================================================================
        Parameter                  Description    Method
        =========================  =============  ========================================================================================
        energy_balance_constraint  :code:`None`   :py:meth:`energy_balance_func <tespy.components.power.bus.PowerBus.energy_balance_func>`
        =========================  =============  ========================================================================================

        Table of parameters

        ===========  =============  ============  ============
        Parameter    Description    Quantity      Method
        ===========  =============  ============  ============
        num_in       :code:`None`   :code:`None`  :code:`None`
        num_out      :code:`None`   :code:`None`  :code:`None`
        ===========  =============  ============  ============




        .. rubric:: Generator

        Class documentation and example: :py:class:`Generator <tespy.components.power.generator.Generator>`

        Table of constraints

        =====================  =============  =======================================================================================================================
        Parameter              Description    Method
        =====================  =============  =======================================================================================================================
        mass_flow_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        =====================  =============  =======================================================================================================================

        Table of parameters

        ===========  =============  ==========  =========================================================================================
        Parameter    Description    Quantity    Method
        ===========  =============  ==========  =========================================================================================
        eta          :code:`None`   efficiency  :py:meth:`eta_func <tespy.components.power.generator.Generator.eta_func>`
        delta_power  :code:`None`   power       :py:meth:`delta_power_func <tespy.components.power.generator.Generator.delta_power_func>`
        ===========  =============  ==========  =========================================================================================

        Table of characteristic lines and maps

        ===========  =============  ===================================================================================
        Parameter    Description    Method
        ===========  =============  ===================================================================================
        eta_char     :code:`None`   :py:meth:`eta_char_func <tespy.components.power.generator.Generator.eta_char_func>`
        ===========  =============  ===================================================================================




        .. rubric:: Motor

        Class documentation and example: :py:class:`Motor <tespy.components.power.motor.Motor>`

        Table of constraints

        =====================  =============  =======================================================================================================================
        Parameter              Description    Method
        =====================  =============  =======================================================================================================================
        mass_flow_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        =====================  =============  =======================================================================================================================

        Table of parameters

        ===========  =============  ==========  =================================================================================
        Parameter    Description    Quantity    Method
        ===========  =============  ==========  =================================================================================
        eta          :code:`None`   efficiency  :py:meth:`eta_func <tespy.components.power.motor.Motor.eta_func>`
        delta_power  :code:`None`   power       :py:meth:`delta_power_func <tespy.components.power.motor.Motor.delta_power_func>`
        ===========  =============  ==========  =================================================================================

        Table of characteristic lines and maps

        ===========  =============  ===========================================================================
        Parameter    Description    Method
        ===========  =============  ===========================================================================
        eta_char     :code:`None`   :py:meth:`eta_char_func <tespy.components.power.motor.Motor.eta_char_func>`
        ===========  =============  ===========================================================================




        .. rubric:: PowerSink

        Class documentation and example: :py:class:`PowerSink <tespy.components.power.sink.PowerSink>`




        .. rubric:: PowerSource

        Class documentation and example: :py:class:`PowerSource <tespy.components.power.source.PowerSource>`



    .. tab-item:: reactors

        .. rubric:: FuelCell

        Class documentation and example: :py:class:`FuelCell <tespy.components.reactors.fuel_cell.FuelCell>`

        Table of constraints

        =============================  =============  ===============================================================================================================================
        Parameter                      Description    Method
        =============================  =============  ===============================================================================================================================
        mass_flow_constraints          :code:`None`   :py:meth:`reactor_mass_flow_func <tespy.components.reactors.fuel_cell.FuelCell.reactor_mass_flow_func>`
        cooling_mass_flow_constraints  :code:`None`   :py:meth:`cooling_mass_flow_structure_matrix <tespy.components.reactors.fuel_cell.FuelCell.cooling_mass_flow_structure_matrix>`
        cooling_fluid_constraints      :code:`None`   :py:meth:`cooling_fluid_structure_matrix <tespy.components.reactors.fuel_cell.FuelCell.cooling_fluid_structure_matrix>`
        energy_balance_constraints     :code:`None`   :py:meth:`energy_balance_func <tespy.components.reactors.fuel_cell.FuelCell.energy_balance_func>`
        reactor_pressure_constraints   :code:`None`   :py:meth:`reactor_pressure_structure_matrix <tespy.components.reactors.fuel_cell.FuelCell.reactor_pressure_structure_matrix>`
        =============================  =============  ===============================================================================================================================

        Table of parameters

        ===========  =============  ===============  ===================================================================================================
        Parameter    Description    Quantity         Method
        ===========  =============  ===============  ===================================================================================================
        P            :code:`None`   power            :code:`None`
        Q            :code:`None`   heat             :py:meth:`heat_func <tespy.components.reactors.fuel_cell.FuelCell.heat_func>`
        pr           :code:`None`   ratio            :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure         :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta         :code:`None`   :code:`None`     :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        eta          :code:`None`   efficiency       :py:meth:`eta_func <tespy.components.reactors.fuel_cell.FuelCell.eta_func>`
        e            :code:`None`   specific_energy  :py:meth:`specific_energy_func <tespy.components.reactors.fuel_cell.FuelCell.specific_energy_func>`
        ===========  =============  ===============  ===================================================================================================




        .. rubric:: WaterElectrolyzer

        Class documentation and example: :py:class:`WaterElectrolyzer <tespy.components.reactors.water_electrolyzer.WaterElectrolyzer>`

        Table of constraints

        =============================  =============  =================================================================================================================================================
        Parameter                      Description    Method
        =============================  =============  =================================================================================================================================================
        mass_flow_constraints          :code:`None`   :py:meth:`reactor_mass_flow_func <tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.reactor_mass_flow_func>`
        cooling_mass_flow_constraints  :code:`None`   :py:meth:`cooling_mass_flow_structure_matrix <tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.cooling_mass_flow_structure_matrix>`
        cooling_fluid_constraints      :code:`None`   :py:meth:`cooling_fluid_structure_matrix <tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.cooling_fluid_structure_matrix>`
        energy_balance_constraints     :code:`None`   :py:meth:`energy_balance_func <tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.energy_balance_func>`
        reactor_pressure_constraints   :code:`None`   :py:meth:`reactor_pressure_structure_matrix <tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.reactor_pressure_structure_matrix>`
        gas_temperature_constraints    :code:`None`   :py:meth:`gas_temperature_func <tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.gas_temperature_func>`
        =============================  =============  =================================================================================================================================================

        Table of parameters

        ===========  =============  ===============  =====================================================================================================================
        Parameter    Description    Quantity         Method
        ===========  =============  ===============  =====================================================================================================================
        P            :code:`None`   power            :code:`None`
        Q            :code:`None`   heat             :py:meth:`heat_func <tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.heat_func>`
        pr           :code:`None`   ratio            :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure         :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta         :code:`None`   :code:`None`     :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        eta          :code:`None`   efficiency       :py:meth:`eta_func <tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.eta_func>`
        e            :code:`None`   specific_energy  :py:meth:`specific_energy_func <tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.specific_energy_func>`
        ===========  =============  ===============  =====================================================================================================================

        Table of characteristic lines and maps

        ===========  =============  =======================================================================================================
        Parameter    Description    Method
        ===========  =============  =======================================================================================================
        eta_char     :code:`None`   :py:meth:`eta_char_func <tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.eta_char_func>`
        ===========  =============  =======================================================================================================



    .. tab-item:: turbomachinery

        .. rubric:: Turbomachine

        Class documentation and example: :py:class:`Turbomachine <tespy.components.turbomachinery.base.Turbomachine>`

        Table of constraints

        =====================  =============  =======================================================================================================================
        Parameter              Description    Method
        =====================  =============  =======================================================================================================================
        mass_flow_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        =====================  =============  =======================================================================================================================

        Table of parameters

        ===========  =============  ==========  ======================================================================================================
        Parameter    Description    Quantity    Method
        ===========  =============  ==========  ======================================================================================================
        P            :code:`None`   power       :py:meth:`energy_balance_func <tespy.components.turbomachinery.base.Turbomachine.energy_balance_func>`
        pr           :code:`None`   ratio       :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure    :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        ===========  =============  ==========  ======================================================================================================




        .. rubric:: Compressor

        Class documentation and example: :py:class:`Compressor <tespy.components.turbomachinery.compressor.Compressor>`

        Table of constraints

        =====================  =============  =======================================================================================================================
        Parameter              Description    Method
        =====================  =============  =======================================================================================================================
        mass_flow_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        =====================  =============  =======================================================================================================================

        Table of parameters

        ===========  =============  ==========  ======================================================================================================
        Parameter    Description    Quantity    Method
        ===========  =============  ==========  ======================================================================================================
        P            :code:`None`   power       :py:meth:`energy_balance_func <tespy.components.turbomachinery.base.Turbomachine.energy_balance_func>`
        pr           :code:`None`   ratio       :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure    :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        eta_s        :code:`None`   efficiency  :py:meth:`eta_s_func <tespy.components.turbomachinery.compressor.Compressor.eta_s_func>`
        igva         :code:`None`   angle       :code:`None`
        ===========  =============  ==========  ======================================================================================================

        Table of parameter groups

        ====================  =============  ====================================  ==========================================================================================================
        Parameter             Description    Required parameters                   Method
        ====================  =============  ====================================  ==========================================================================================================
        char_map_eta_s_group  :code:`None`   :code:`char_map_eta_s`, :code:`igva`  :py:meth:`char_map_eta_s_func <tespy.components.turbomachinery.compressor.Compressor.char_map_eta_s_func>`
        char_map_pr_group     :code:`None`   :code:`char_map_pr`, :code:`igva`     :py:meth:`char_map_pr_func <tespy.components.turbomachinery.compressor.Compressor.char_map_pr_func>`
        ====================  =============  ====================================  ==========================================================================================================

        Table of characteristic lines and maps

        ==============  =============  ==================================================================================================
        Parameter       Description    Method
        ==============  =============  ==================================================================================================
        eta_s_char      :code:`None`   :py:meth:`eta_s_char_func <tespy.components.turbomachinery.compressor.Compressor.eta_s_char_func>`
        char_map_eta_s  :code:`None`   :code:`None`
        char_map_pr     :code:`None`   :code:`None`
        ==============  =============  ==================================================================================================




        .. rubric:: Pump

        Class documentation and example: :py:class:`Pump <tespy.components.turbomachinery.pump.Pump>`

        Table of constraints

        =====================  =============  =======================================================================================================================
        Parameter              Description    Method
        =====================  =============  =======================================================================================================================
        mass_flow_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        =====================  =============  =======================================================================================================================

        Table of parameters

        ===========  =============  ==========  ======================================================================================================
        Parameter    Description    Quantity    Method
        ===========  =============  ==========  ======================================================================================================
        P            :code:`None`   power       :py:meth:`energy_balance_func <tespy.components.turbomachinery.base.Turbomachine.energy_balance_func>`
        pr           :code:`None`   ratio       :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure    :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        eta_s        :code:`None`   efficiency  :py:meth:`eta_s_func <tespy.components.turbomachinery.pump.Pump.eta_s_func>`
        ===========  =============  ==========  ======================================================================================================

        Table of characteristic lines and maps

        ===========  =============  ======================================================================================
        Parameter    Description    Method
        ===========  =============  ======================================================================================
        eta_s_char   :code:`None`   :py:meth:`eta_s_char_func <tespy.components.turbomachinery.pump.Pump.eta_s_char_func>`
        flow_char    :code:`None`   :py:meth:`flow_char_func <tespy.components.turbomachinery.pump.Pump.flow_char_func>`
        ===========  =============  ======================================================================================




        .. rubric:: Turbine

        Class documentation and example: :py:class:`Turbine <tespy.components.turbomachinery.turbine.Turbine>`

        Table of constraints

        =====================  =============  =======================================================================================================================
        Parameter              Description    Method
        =====================  =============  =======================================================================================================================
        mass_flow_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        =====================  =============  =======================================================================================================================

        Table of parameters

        ===========  =============  ============  ======================================================================================================
        Parameter    Description    Quantity      Method
        ===========  =============  ============  ======================================================================================================
        P            :code:`None`   power         :py:meth:`energy_balance_func <tespy.components.turbomachinery.base.Turbomachine.energy_balance_func>`
        pr           :code:`None`   ratio         :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure      :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        eta_s        :code:`None`   efficiency    :py:meth:`eta_s_func <tespy.components.turbomachinery.turbine.Turbine.eta_s_func>`
        cone         :code:`None`   :code:`None`  :py:meth:`cone_func <tespy.components.turbomachinery.turbine.Turbine.cone_func>`
        ===========  =============  ============  ======================================================================================================

        Table of characteristic lines and maps

        ===========  =============  ============================================================================================
        Parameter    Description    Method
        ===========  =============  ============================================================================================
        eta_s_char   :code:`None`   :py:meth:`eta_s_char_func <tespy.components.turbomachinery.turbine.Turbine.eta_s_char_func>`
        ===========  =============  ============================================================================================




        .. rubric:: SteamTurbine

        Class documentation and example: :py:class:`SteamTurbine <tespy.components.turbomachinery.steam_turbine.SteamTurbine>`

        Table of constraints

        =====================  =============  =======================================================================================================================
        Parameter              Description    Method
        =====================  =============  =======================================================================================================================
        mass_flow_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        =====================  =============  =======================================================================================================================

        Table of parameters

        ===========  =============  ============  ======================================================================================================
        Parameter    Description    Quantity      Method
        ===========  =============  ============  ======================================================================================================
        P            :code:`None`   power         :py:meth:`energy_balance_func <tespy.components.turbomachinery.base.Turbomachine.energy_balance_func>`
        pr           :code:`None`   ratio         :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure      :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        eta_s        :code:`None`   efficiency    :py:meth:`eta_s_func <tespy.components.turbomachinery.turbine.Turbine.eta_s_func>`
        cone         :code:`None`   :code:`None`  :py:meth:`cone_func <tespy.components.turbomachinery.turbine.Turbine.cone_func>`
        alpha        :code:`None`   ratio         :code:`None`
        eta_s_dry    :code:`None`   efficiency    :code:`None`
        ===========  =============  ============  ======================================================================================================

        Table of parameter groups

        ===============  =============  ================================  =====================================================================================================
        Parameter        Description    Required parameters               Method
        ===============  =============  ================================  =====================================================================================================
        eta_s_dry_group  :code:`None`   :code:`alpha`, :code:`eta_s_dry`  :py:meth:`eta_s_wet_func <tespy.components.turbomachinery.steam_turbine.SteamTurbine.eta_s_wet_func>`
        ===============  =============  ================================  =====================================================================================================

        Table of characteristic lines and maps

        ===========  =============  ============================================================================================
        Parameter    Description    Method
        ===========  =============  ============================================================================================
        eta_s_char   :code:`None`   :py:meth:`eta_s_char_func <tespy.components.turbomachinery.turbine.Turbine.eta_s_char_func>`
        ===========  =============  ============================================================================================




        .. rubric:: TurboCompressor

        Class documentation and example: :py:class:`TurboCompressor <tespy.components.turbomachinery.turbocompressor.TurboCompressor>`

        Table of constraints

        =====================  =============  =======================================================================================================================
        Parameter              Description    Method
        =====================  =============  =======================================================================================================================
        mass_flow_constraints  :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        fluid_constraints      :code:`None`   :py:meth:`variable_equality_structure_matrix <tespy.components.component.Component.variable_equality_structure_matrix>`
        =====================  =============  =======================================================================================================================

        Table of parameters

        ===========  =============  ==========  ======================================================================================================
        Parameter    Description    Quantity    Method
        ===========  =============  ==========  ======================================================================================================
        P            :code:`None`   power       :py:meth:`energy_balance_func <tespy.components.turbomachinery.base.Turbomachine.energy_balance_func>`
        pr           :code:`None`   ratio       :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure    :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        eta_s        :code:`None`   efficiency  :py:meth:`eta_s_func <tespy.components.turbomachinery.compressor.Compressor.eta_s_func>`
        igva         :code:`None`   angle       :code:`None`
        ===========  =============  ==========  ======================================================================================================

        Table of parameter groups

        ====================  =============  ====================================  ====================================================================================================================
        Parameter             Description    Required parameters                   Method
        ====================  =============  ====================================  ====================================================================================================================
        char_map_eta_s_group  :code:`None`   :code:`char_map_eta_s`, :code:`igva`  :py:meth:`char_map_eta_s_func <tespy.components.turbomachinery.turbocompressor.TurboCompressor.char_map_eta_s_func>`
        char_map_pr_group     :code:`None`   :code:`char_map_pr`, :code:`igva`     :py:meth:`char_map_pr_func <tespy.components.turbomachinery.turbocompressor.TurboCompressor.char_map_pr_func>`
        ====================  =============  ====================================  ====================================================================================================================

        Table of characteristic lines and maps

        ==============  =============  ============
        Parameter       Description    Method
        ==============  =============  ============
        char_map_eta_s  :code:`None`   :code:`None`
        char_map_pr     :code:`None`   :code:`None`
        ==============  =============  ============
