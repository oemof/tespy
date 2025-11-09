.. tab-set::

    .. tab-item:: components

        .. rubric:: Component
        
        Class documentation and example: :py:class:`Component <tespy.components.component.Component>`
        
        

    .. tab-item:: basics

        .. rubric:: CycleCloser
        
        Class documentation and example: :py:class:`CycleCloser <tespy.components.basics.cycle_closer.CycleCloser>`
        
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
        
        


        .. rubric:: SubsystemInterface
        
        Class documentation and example: :py:class:`SubsystemInterface <tespy.components.basics.subsystem_interface.SubsystemInterface>`
        
        ===========  =============  ============  ============
        Parameter    Description    Quantity      Method
        ===========  =============  ============  ============
        num_inter    :code:`None`   :code:`None`  :code:`None`
        ===========  =============  ============  ============
        
        

    .. tab-item:: combustion

        .. rubric:: CombustionChamber
        
        Class documentation and example: :py:class:`CombustionChamber <tespy.components.combustion.base.CombustionChamber>`
        
        ===========  =============  ==========  =======================================================================================
        Parameter    Description    Quantity    Method
        ===========  =============  ==========  =======================================================================================
        lamb         :code:`None`   ratio       :py:meth:`lambda_func <tespy.components.combustion.base.CombustionChamber.lambda_func>`
        ti           :code:`None`   heat        :py:meth:`ti_func <tespy.components.combustion.base.CombustionChamber.ti_func>`
        ===========  =============  ==========  =======================================================================================
        
        


        .. rubric:: DiabaticCombustionChamber
        
        Class documentation and example: :py:class:`DiabaticCombustionChamber <tespy.components.combustion.diabatic.DiabaticCombustionChamber>`
        
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
        
        ===========  =============  ============  =========================================================================================
        Parameter    Description    Quantity      Method
        ===========  =============  ============  =========================================================================================
        lamb         :code:`None`   ratio         :py:meth:`lambda_func <tespy.components.combustion.base.CombustionChamber.lambda_func>`
        ti           :code:`None`   heat          :py:meth:`ti_func <tespy.components.combustion.base.CombustionChamber.ti_func>`
        Q1           :code:`None`   heat          :py:meth:`Q1_func <tespy.components.combustion.engine.CombustionEngine.Q1_func>`
        Q2           :code:`None`   heat          :py:meth:`Q2_func <tespy.components.combustion.engine.CombustionEngine.Q2_func>`
        pr1          :code:`None`   ratio         :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        pr2          :code:`None`   ratio         :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp1          :code:`None`   pressure      :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        dp2          :code:`None`   pressure      :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta1        :code:`None`   :code:`None`  :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        zeta2        :code:`None`   :code:`None`  :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        P            :code:`None`   power         :code:`None`
        Qloss        :code:`None`   heat          :code:`None`
        eta_mech     :code:`None`   :code:`None`  :code:`None`
        T_v_inner    :code:`None`   :code:`None`  :code:`None`
        ===========  =============  ============  =========================================================================================
        
        

    .. tab-item:: displacementmachinery

        .. rubric:: DisplacementMachine
        
        Class documentation and example: :py:class:`DisplacementMachine <tespy.components.displacementmachinery.base.DisplacementMachine>`
        
        ===========  =============  ==========  ====================================================================================================================
        Parameter    Description    Quantity    Method
        ===========  =============  ==========  ====================================================================================================================
        P            :code:`None`   power       :py:meth:`energy_balance_func <tespy.components.displacementmachinery.base.DisplacementMachine.energy_balance_func>`
        pr           :code:`None`   ratio       :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure    :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        ===========  =============  ==========  ====================================================================================================================
        
        


        .. rubric:: PolynomialCompressor
        
        Class documentation and example: :py:class:`PolynomialCompressor <tespy.components.displacementmachinery.polynomial_compressor.PolynomialCompressor>`
        
        =================  =============  ============  =========================================================================================
        Parameter          Description    Quantity      Method
        =================  =============  ============  =========================================================================================
        pr                 :code:`None`   ratio         :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp                 :code:`None`   pressure      :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        P                  :code:`None`   power         :code:`None`
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
        
        


        .. rubric:: PolynomialCompressorWithCooling
        
        Class documentation and example: :py:class:`PolynomialCompressorWithCooling <tespy.components.displacementmachinery.polynomial_compressor_with_cooling.PolynomialCompressorWithCooling>`
        
        =================  =============  ======================  =========================================================================================
        Parameter          Description    Quantity                Method
        =================  =============  ======================  =========================================================================================
        pr                 :code:`None`   ratio                   :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp                 :code:`None`   pressure                :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        dp_cooling         :code:`None`   pressure                :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        pr_cooling         :code:`None`   ratio                   :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        P                  :code:`None`   power                   :code:`None`
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
        =================  =============  ======================  =========================================================================================
        
        

    .. tab-item:: heat_exchangers

        .. rubric:: HeatExchanger
        
        Class documentation and example: :py:class:`HeatExchanger <tespy.components.heat_exchangers.base.HeatExchanger>`
        
        ===========  =============  =========================  ================================================================================================================
        Parameter    Description    Quantity                   Method
        ===========  =============  =========================  ================================================================================================================
        Q            :code:`None`   heat                       :py:meth:`energy_balance_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func>`
        kA           :code:`None`   heat_transfer_coefficient  :py:meth:`kA_func <tespy.components.heat_exchangers.base.HeatExchanger.kA_func>`
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
        td_log       :code:`None`   temperature_difference     :code:`None`
        ===========  =============  =========================  ================================================================================================================
        
        


        .. rubric:: Condenser
        
        Class documentation and example: :py:class:`Condenser <tespy.components.heat_exchangers.condenser.Condenser>`
        
        ===========  =============  =========================  ================================================================================================================
        Parameter    Description    Quantity                   Method
        ===========  =============  =========================  ================================================================================================================
        Q            :code:`None`   heat                       :py:meth:`energy_balance_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func>`
        kA           :code:`None`   heat_transfer_coefficient  :py:meth:`kA_func <tespy.components.heat_exchangers.base.HeatExchanger.kA_func>`
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
        td_log       :code:`None`   temperature_difference     :code:`None`
        ===========  =============  =========================  ================================================================================================================
        
        


        .. rubric:: Desuperheater
        
        Class documentation and example: :py:class:`Desuperheater <tespy.components.heat_exchangers.desuperheater.Desuperheater>`
        
        ===========  =============  =========================  ================================================================================================================
        Parameter    Description    Quantity                   Method
        ===========  =============  =========================  ================================================================================================================
        Q            :code:`None`   heat                       :py:meth:`energy_balance_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func>`
        kA           :code:`None`   heat_transfer_coefficient  :py:meth:`kA_func <tespy.components.heat_exchangers.base.HeatExchanger.kA_func>`
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
        td_log       :code:`None`   temperature_difference     :code:`None`
        ===========  =============  =========================  ================================================================================================================
        
        


        .. rubric:: SectionedHeatExchanger
        
        Class documentation and example: :py:class:`SectionedHeatExchanger <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger>`
        
        =================  =============  =========================  ================================================================================================================
        Parameter          Description    Quantity                   Method
        =================  =============  =========================  ================================================================================================================
        Q                  :code:`None`   heat                       :py:meth:`energy_balance_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func>`
        kA                 :code:`None`   heat_transfer_coefficient  :py:meth:`kA_func <tespy.components.heat_exchangers.base.HeatExchanger.kA_func>`
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
        td_pinch           :code:`None`   temperature_difference     :py:meth:`td_pinch_func <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.td_pinch_func>`
        td_log             :code:`None`   temperature_difference     :code:`None`
        num_sections       :code:`None`   :code:`None`               :code:`None`
        refrigerant_index  :code:`None`   :code:`None`               :code:`None`
        re_exp_r           :code:`None`   :code:`None`               :code:`None`
        re_exp_sf          :code:`None`   :code:`None`               :code:`None`
        alpha_ratio        :code:`None`   ratio                      :code:`None`
        area_ratio         :code:`None`   ratio                      :code:`None`
        =================  =============  =========================  ================================================================================================================
        
        


        .. rubric:: MovingBoundaryHeatExchanger
        
        Class documentation and example: :py:class:`MovingBoundaryHeatExchanger <tespy.components.heat_exchangers.movingboundary.MovingBoundaryHeatExchanger>`
        
        =================  =============  =========================  ================================================================================================================
        Parameter          Description    Quantity                   Method
        =================  =============  =========================  ================================================================================================================
        Q                  :code:`None`   heat                       :py:meth:`energy_balance_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func>`
        kA                 :code:`None`   heat_transfer_coefficient  :py:meth:`kA_func <tespy.components.heat_exchangers.base.HeatExchanger.kA_func>`
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
        td_pinch           :code:`None`   temperature_difference     :py:meth:`td_pinch_func <tespy.components.heat_exchangers.sectioned.SectionedHeatExchanger.td_pinch_func>`
        td_log             :code:`None`   temperature_difference     :code:`None`
        refrigerant_index  :code:`None`   :code:`None`               :code:`None`
        re_exp_r           :code:`None`   :code:`None`               :code:`None`
        re_exp_sf          :code:`None`   :code:`None`               :code:`None`
        alpha_ratio        :code:`None`   ratio                      :code:`None`
        area_ratio         :code:`None`   ratio                      :code:`None`
        =================  =============  =========================  ================================================================================================================
        
        


        .. rubric:: SimpleHeatExchanger
        
        Class documentation and example: :py:class:`SimpleHeatExchanger <tespy.components.heat_exchangers.simple.SimpleHeatExchanger>`
        
        ========================  =============  =========================  ================================================================================================================
        Parameter                 Description    Quantity                   Method
        ========================  =============  =========================  ================================================================================================================
        Q                         :code:`None`   heat                       :py:meth:`energy_balance_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_balance_func>`
        pr                        :code:`None`   ratio                      :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp                        :code:`None`   pressure                   :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta                      :code:`None`   :code:`None`               :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        power_connector_location  :code:`None`   :code:`None`               :code:`None`
        D                         :code:`None`   length                     :code:`None`
        L                         :code:`None`   length                     :code:`None`
        ks                        :code:`None`   length                     :code:`None`
        ks_HW                     :code:`None`   :code:`None`               :code:`None`
        kA                        :code:`None`   heat_transfer_coefficient  :code:`None`
        Tamb                      :code:`None`   temperature                :code:`None`
        dissipative               :code:`None`   :code:`None`               :code:`None`
        ========================  =============  =========================  ================================================================================================================
        
        


        .. rubric:: ParabolicTrough
        
        Class documentation and example: :py:class:`ParabolicTrough <tespy.components.heat_exchangers.parabolic_trough.ParabolicTrough>`
        
        ========================  =============  ============  ================================================================================================================
        Parameter                 Description    Quantity      Method
        ========================  =============  ============  ================================================================================================================
        Q                         :code:`None`   heat          :py:meth:`energy_balance_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_balance_func>`
        pr                        :code:`None`   ratio         :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp                        :code:`None`   pressure      :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta                      :code:`None`   :code:`None`  :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        power_connector_location  :code:`None`   :code:`None`  :code:`None`
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
        
        


        .. rubric:: ParallelFlowHeatExchanger
        
        Class documentation and example: :py:class:`ParallelFlowHeatExchanger <tespy.components.heat_exchangers.parallel.ParallelFlowHeatExchanger>`
        
        ===========  =============  =========================  ================================================================================================================
        Parameter    Description    Quantity                   Method
        ===========  =============  =========================  ================================================================================================================
        Q            :code:`None`   heat                       :py:meth:`energy_balance_hot_func <tespy.components.heat_exchangers.base.HeatExchanger.energy_balance_hot_func>`
        kA           :code:`None`   heat_transfer_coefficient  :py:meth:`kA_func <tespy.components.heat_exchangers.base.HeatExchanger.kA_func>`
        ttd_u        :code:`None`   temperature_difference     :py:meth:`ttd_u_func <tespy.components.heat_exchangers.parallel.ParallelFlowHeatExchanger.ttd_u_func>`
        ttd_l        :code:`None`   temperature_difference     :py:meth:`ttd_l_func <tespy.components.heat_exchangers.parallel.ParallelFlowHeatExchanger.ttd_l_func>`
        pr1          :code:`None`   ratio                      :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        pr2          :code:`None`   ratio                      :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp1          :code:`None`   pressure                   :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        dp2          :code:`None`   pressure                   :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta1        :code:`None`   :code:`None`               :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        zeta2        :code:`None`   :code:`None`               :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        td_log       :code:`None`   temperature_difference     :code:`None`
        ===========  =============  =========================  ================================================================================================================
        
        


        .. rubric:: SolarCollector
        
        Class documentation and example: :py:class:`SolarCollector <tespy.components.heat_exchangers.solar_collector.SolarCollector>`
        
        ========================  =============  ============  ================================================================================================================
        Parameter                 Description    Quantity      Method
        ========================  =============  ============  ================================================================================================================
        Q                         :code:`None`   heat          :py:meth:`energy_balance_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_balance_func>`
        pr                        :code:`None`   ratio         :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp                        :code:`None`   pressure      :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta                      :code:`None`   :code:`None`  :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        power_connector_location  :code:`None`   :code:`None`  :code:`None`
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
        
        

    .. tab-item:: nodes

        .. rubric:: NodeBase
        
        Class documentation and example: :py:class:`NodeBase <tespy.components.nodes.base.NodeBase>`
        
        


        .. rubric:: DropletSeparator
        
        Class documentation and example: :py:class:`DropletSeparator <tespy.components.nodes.droplet_separator.DropletSeparator>`
        
        


        .. rubric:: Drum
        
        Class documentation and example: :py:class:`Drum <tespy.components.nodes.drum.Drum>`
        
        


        .. rubric:: Merge
        
        Class documentation and example: :py:class:`Merge <tespy.components.nodes.merge.Merge>`
        
        ===========  =============  ============  ============
        Parameter    Description    Quantity      Method
        ===========  =============  ============  ============
        num_in       :code:`None`   :code:`None`  :code:`None`
        ===========  =============  ============  ============
        
        


        .. rubric:: Splitter
        
        Class documentation and example: :py:class:`Splitter <tespy.components.nodes.splitter.Splitter>`
        
        ===========  =============  ============  ============
        Parameter    Description    Quantity      Method
        ===========  =============  ============  ============
        num_out      :code:`None`   :code:`None`  :code:`None`
        ===========  =============  ============  ============
        
        


        .. rubric:: Node
        
        Class documentation and example: :py:class:`Node <tespy.components.nodes.node.Node>`
        
        ===========  =============  ============  ============
        Parameter    Description    Quantity      Method
        ===========  =============  ============  ============
        num_out      :code:`None`   :code:`None`  :code:`None`
        num_in       :code:`None`   :code:`None`  :code:`None`
        ===========  =============  ============  ============
        
        


        .. rubric:: Separator
        
        Class documentation and example: :py:class:`Separator <tespy.components.nodes.separator.Separator>`
        
        ===========  =============  ============  ============
        Parameter    Description    Quantity      Method
        ===========  =============  ============  ============
        num_out      :code:`None`   :code:`None`  :code:`None`
        ===========  =============  ============  ============
        
        

    .. tab-item:: piping

        .. rubric:: Pipe
        
        Class documentation and example: :py:class:`Pipe <tespy.components.piping.pipe.Pipe>`
        
        ========================  =============  =========================  ================================================================================================================
        Parameter                 Description    Quantity                   Method
        ========================  =============  =========================  ================================================================================================================
        Q                         :code:`None`   heat                       :py:meth:`energy_balance_func <tespy.components.heat_exchangers.simple.SimpleHeatExchanger.energy_balance_func>`
        pr                        :code:`None`   ratio                      :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp                        :code:`None`   pressure                   :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta                      :code:`None`   :code:`None`               :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        power_connector_location  :code:`None`   :code:`None`               :code:`None`
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
        
        


        .. rubric:: Valve
        
        Class documentation and example: :py:class:`Valve <tespy.components.piping.valve.Valve>`
        
        ===========  =============  ============  =========================================================================================
        Parameter    Description    Quantity      Method
        ===========  =============  ============  =========================================================================================
        pr           :code:`None`   ratio         :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure      :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta         :code:`None`   :code:`None`  :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        ===========  =============  ============  =========================================================================================
        
        

    .. tab-item:: power

        .. rubric:: PowerBus
        
        Class documentation and example: :py:class:`PowerBus <tespy.components.power.bus.PowerBus>`
        
        ===========  =============  ============  ============
        Parameter    Description    Quantity      Method
        ===========  =============  ============  ============
        num_in       :code:`None`   :code:`None`  :code:`None`
        num_out      :code:`None`   :code:`None`  :code:`None`
        ===========  =============  ============  ============
        
        


        .. rubric:: Generator
        
        Class documentation and example: :py:class:`Generator <tespy.components.power.generator.Generator>`
        
        ===========  =============  ==========  =========================================================================================
        Parameter    Description    Quantity    Method
        ===========  =============  ==========  =========================================================================================
        eta          :code:`None`   efficiency  :py:meth:`eta_func <tespy.components.power.generator.Generator.eta_func>`
        delta_power  :code:`None`   power       :py:meth:`delta_power_func <tespy.components.power.generator.Generator.delta_power_func>`
        ===========  =============  ==========  =========================================================================================
        
        


        .. rubric:: Motor
        
        Class documentation and example: :py:class:`Motor <tespy.components.power.motor.Motor>`
        
        ===========  =============  ==========  =================================================================================
        Parameter    Description    Quantity    Method
        ===========  =============  ==========  =================================================================================
        eta          :code:`None`   efficiency  :py:meth:`eta_func <tespy.components.power.motor.Motor.eta_func>`
        delta_power  :code:`None`   power       :py:meth:`delta_power_func <tespy.components.power.motor.Motor.delta_power_func>`
        ===========  =============  ==========  =================================================================================
        
        


        .. rubric:: PowerSink
        
        Class documentation and example: :py:class:`PowerSink <tespy.components.power.sink.PowerSink>`
        
        


        .. rubric:: PowerSource
        
        Class documentation and example: :py:class:`PowerSource <tespy.components.power.source.PowerSource>`
        
        

    .. tab-item:: reactors

        .. rubric:: FuelCell
        
        Class documentation and example: :py:class:`FuelCell <tespy.components.reactors.fuel_cell.FuelCell>`
        
        ===========  =============  ===============  ===================================================================================================
        Parameter    Description    Quantity         Method
        ===========  =============  ===============  ===================================================================================================
        Q            :code:`None`   heat             :py:meth:`heat_func <tespy.components.reactors.fuel_cell.FuelCell.heat_func>`
        pr           :code:`None`   ratio            :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure         :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta         :code:`None`   :code:`None`     :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        eta          :code:`None`   efficiency       :py:meth:`eta_func <tespy.components.reactors.fuel_cell.FuelCell.eta_func>`
        e            :code:`None`   specific_energy  :py:meth:`specific_energy_func <tespy.components.reactors.fuel_cell.FuelCell.specific_energy_func>`
        P            :code:`None`   power            :code:`None`
        ===========  =============  ===============  ===================================================================================================
        
        


        .. rubric:: WaterElectrolyzer
        
        Class documentation and example: :py:class:`WaterElectrolyzer <tespy.components.reactors.water_electrolyzer.WaterElectrolyzer>`
        
        ===========  =============  ===============  =====================================================================================================================
        Parameter    Description    Quantity         Method
        ===========  =============  ===============  =====================================================================================================================
        Q            :code:`None`   heat             :py:meth:`heat_func <tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.heat_func>`
        pr           :code:`None`   ratio            :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure         :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        zeta         :code:`None`   :code:`None`     :py:meth:`zeta_func <tespy.components.component.Component.zeta_func>`
        eta          :code:`None`   efficiency       :py:meth:`eta_func <tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.eta_func>`
        e            :code:`None`   specific_energy  :py:meth:`specific_energy_func <tespy.components.reactors.water_electrolyzer.WaterElectrolyzer.specific_energy_func>`
        P            :code:`None`   power            :code:`None`
        ===========  =============  ===============  =====================================================================================================================
        
        

    .. tab-item:: turbomachinery

        .. rubric:: Turbomachine
        
        Class documentation and example: :py:class:`Turbomachine <tespy.components.turbomachinery.base.Turbomachine>`
        
        ===========  =============  ==========  ======================================================================================================
        Parameter    Description    Quantity    Method
        ===========  =============  ==========  ======================================================================================================
        P            :code:`None`   power       :py:meth:`energy_balance_func <tespy.components.turbomachinery.base.Turbomachine.energy_balance_func>`
        pr           :code:`None`   ratio       :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure    :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        ===========  =============  ==========  ======================================================================================================
        
        


        .. rubric:: Compressor
        
        Class documentation and example: :py:class:`Compressor <tespy.components.turbomachinery.compressor.Compressor>`
        
        ===========  =============  ==========  ======================================================================================================
        Parameter    Description    Quantity    Method
        ===========  =============  ==========  ======================================================================================================
        P            :code:`None`   power       :py:meth:`energy_balance_func <tespy.components.turbomachinery.base.Turbomachine.energy_balance_func>`
        pr           :code:`None`   ratio       :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure    :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        eta_s        :code:`None`   efficiency  :py:meth:`eta_s_func <tespy.components.turbomachinery.compressor.Compressor.eta_s_func>`
        igva         :code:`None`   angle       :code:`None`
        ===========  =============  ==========  ======================================================================================================
        
        


        .. rubric:: Pump
        
        Class documentation and example: :py:class:`Pump <tespy.components.turbomachinery.pump.Pump>`
        
        ===========  =============  ==========  ======================================================================================================
        Parameter    Description    Quantity    Method
        ===========  =============  ==========  ======================================================================================================
        P            :code:`None`   power       :py:meth:`energy_balance_func <tespy.components.turbomachinery.base.Turbomachine.energy_balance_func>`
        pr           :code:`None`   ratio       :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure    :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        eta_s        :code:`None`   efficiency  :py:meth:`eta_s_func <tespy.components.turbomachinery.pump.Pump.eta_s_func>`
        ===========  =============  ==========  ======================================================================================================
        
        


        .. rubric:: Turbine
        
        Class documentation and example: :py:class:`Turbine <tespy.components.turbomachinery.turbine.Turbine>`
        
        ===========  =============  ============  ======================================================================================================
        Parameter    Description    Quantity      Method
        ===========  =============  ============  ======================================================================================================
        P            :code:`None`   power         :py:meth:`energy_balance_func <tespy.components.turbomachinery.base.Turbomachine.energy_balance_func>`
        pr           :code:`None`   ratio         :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure      :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        eta_s        :code:`None`   efficiency    :py:meth:`eta_s_func <tespy.components.turbomachinery.turbine.Turbine.eta_s_func>`
        cone         :code:`None`   :code:`None`  :py:meth:`cone_func <tespy.components.turbomachinery.turbine.Turbine.cone_func>`
        ===========  =============  ============  ======================================================================================================
        
        


        .. rubric:: SteamTurbine
        
        Class documentation and example: :py:class:`SteamTurbine <tespy.components.turbomachinery.steam_turbine.SteamTurbine>`
        
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
        
        


        .. rubric:: TurboCompressor
        
        Class documentation and example: :py:class:`TurboCompressor <tespy.components.turbomachinery.turbocompressor.TurboCompressor>`
        
        ===========  =============  ==========  ======================================================================================================
        Parameter    Description    Quantity    Method
        ===========  =============  ==========  ======================================================================================================
        P            :code:`None`   power       :py:meth:`energy_balance_func <tespy.components.turbomachinery.base.Turbomachine.energy_balance_func>`
        pr           :code:`None`   ratio       :py:meth:`pr_structure_matrix <tespy.components.component.Component.pr_structure_matrix>`
        dp           :code:`None`   pressure    :py:meth:`dp_structure_matrix <tespy.components.component.Component.dp_structure_matrix>`
        eta_s        :code:`None`   efficiency  :py:meth:`eta_s_func <tespy.components.turbomachinery.compressor.Compressor.eta_s_func>`
        igva         :code:`None`   angle       :code:`None`
        ===========  =============  ==========  ======================================================================================================
        
        

