<?xml version='1.0' encoding='ISO-8859-1' standalone='yes'?>
<tagfile>
  <compound kind="page">
    <name>index</name>
    <title>LSMLIB Documentation</title>
    <filename>index</filename>
  </compound>
  <compound kind="file">
    <name>BoundaryConditionModule.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/parallel/</path>
    <filename>BoundaryConditionModule_8h</filename>
    <includes id="LevelSetMethodToolbox_8h" name="LevelSetMethodToolbox.h" local="yes">LevelSetMethodToolbox.h</includes>
    <namespace>hier</namespace>
    <namespace>LSMLIB</namespace>
    <namespace>SAMRAI</namespace>
    <namespace>tbox</namespace>
  </compound>
  <compound kind="file">
    <name>FieldExtensionAlgorithm.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/parallel/</path>
    <filename>FieldExtensionAlgorithm_8h</filename>
    <includes id="BoundaryConditionModule_8h" name="BoundaryConditionModule.h" local="yes">BoundaryConditionModule.h</includes>
    <includes id="LevelSetMethodToolbox_8h" name="LevelSetMethodToolbox.h" local="yes">LevelSetMethodToolbox.h</includes>
    <namespace>LSMLIB</namespace>
    <namespace>xfer</namespace>
  </compound>
  <compound kind="file">
    <name>FMM_Callback_API.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/fast_marching_method/</path>
    <filename>FMM__Callback__API_8h</filename>
    <includes id="FMM__Core_8h" name="FMM_Core.h" local="yes">FMM_Core.h</includes>
    <member kind="typedef">
      <type>FMM_FieldData</type>
      <name>FMM_FieldData</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initializeFront_CallbackFunc</name>
      <anchor>a1</anchor>
      <arglist>(FMM_CoreData *fmm_core_data, FMM_FieldData *fmm_field_data, int num_dims, int *grid_dims, LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="function">
      <type>LSMLIB_REAL</type>
      <name>updateGridPoint_CallbackFunc</name>
      <anchor>a2</anchor>
      <arglist>(FMM_CoreData *fmm_core_data, FMM_FieldData *fmm_field_data, int *grid_idx, int num_dims, int *grid_dims, LSMLIB_REAL *dx)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>FMM_Core.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/fast_marching_method/</path>
    <filename>FMM__Core_8h</filename>
    <member kind="typedef">
      <type>FMM_CoreData</type>
      <name>FMM_CoreData</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>FMM_FieldData</type>
      <name>FMM_FieldData</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>void(*</type>
      <name>initializeFrontFuncPtr</name>
      <anchor>a2</anchor>
      <arglist>)(FMM_CoreData *fmm_core_data, FMM_FieldData *fmm_field_data, int num_dims, int *grid_dims, LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="typedef">
      <type>LSMLIB_REAL(*</type>
      <name>updateGridPointFuncPtr</name>
      <anchor>a3</anchor>
      <arglist>)(FMM_CoreData *fmm_core_data, FMM_FieldData *fmm_field_data, int *grid_idx, int num_dims, int *grid_dims, LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="enumeration">
      <name>PointStatus</name>
      <anchor>a16</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>KNOWN</name>
      <anchor>a16a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>TRIAL</name>
      <anchor>a16a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>FAR</name>
      <anchor>a16a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>OUTSIDE_DOMAIN</name>
      <anchor>a16a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>FMM_CoreData *</type>
      <name>FMM_Core_createFMM_CoreData</name>
      <anchor>a8</anchor>
      <arglist>(FMM_FieldData *fmm_field_data, int num_dims, int *grid_dims, LSMLIB_REAL *dx, initializeFrontFuncPtr initializeFront, updateGridPointFuncPtr updateGridPoint)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FMM_Core_destroyFMM_CoreData</name>
      <anchor>a9</anchor>
      <arglist>(FMM_CoreData *fmm_core_data)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FMM_Core_initializeFront</name>
      <anchor>a10</anchor>
      <arglist>(FMM_CoreData *fmm_core_data)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FMM_Core_setInitialFrontPoint</name>
      <anchor>a11</anchor>
      <arglist>(FMM_CoreData *fmm_core_data, int *grid_idx, LSMLIB_REAL value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FMM_Core_markPointOutsideDomain</name>
      <anchor>a12</anchor>
      <arglist>(FMM_CoreData *fmm_core_data, int *grid_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FMM_Core_advanceFront</name>
      <anchor>a13</anchor>
      <arglist>(FMM_CoreData *fmm_core_data)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>FMM_Core_moreGridPointsToUpdate</name>
      <anchor>a14</anchor>
      <arglist>(FMM_CoreData *fmm_core_data)</arglist>
    </member>
    <member kind="function">
      <type>int *</type>
      <name>FMM_Core_getGridPointStatusDataArray</name>
      <anchor>a15</anchor>
      <arglist>(FMM_CoreData *fmm_core_data)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>FMM_Heap.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/fast_marching_method/</path>
    <filename>FMM__Heap_8h</filename>
    <class kind="struct">HeapNode</class>
    <member kind="define">
      <type>#define</type>
      <name>FMM_HEAP_MAX_NDIM</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>FMM_Heap</type>
      <name>FMM_Heap</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>HeapNode</type>
      <name>FMM_HeapNode</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>FMM_Heap *</type>
      <name>FMM_Heap_createHeap</name>
      <anchor>a3</anchor>
      <arglist>(int num_dims, int heap_mem_size, LSMLIB_REAL growth_factor)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FMM_Heap_destroyHeap</name>
      <anchor>a4</anchor>
      <arglist>(FMM_Heap *heap)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>FMM_Heap_insertNode</name>
      <anchor>a5</anchor>
      <arglist>(FMM_Heap *heap, int *grid_idx, LSMLIB_REAL value)</arglist>
    </member>
    <member kind="function">
      <type>FMM_HeapNode</type>
      <name>FMM_Heap_extractMin</name>
      <anchor>a6</anchor>
      <arglist>(FMM_Heap *heap, FMM_HeapNode *moved_node, int *moved_handle)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FMM_Heap_updateNode</name>
      <anchor>a7</anchor>
      <arglist>(FMM_Heap *heap, int node_handle, LSMLIB_REAL value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FMM_Heap_clear</name>
      <anchor>a8</anchor>
      <arglist>(FMM_Heap *heap)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>FMM_Heap_isEmpty</name>
      <anchor>a9</anchor>
      <arglist>(FMM_Heap *heap)</arglist>
    </member>
    <member kind="function">
      <type>FMM_HeapNode</type>
      <name>FMM_Heap_getNode</name>
      <anchor>a10</anchor>
      <arglist>(FMM_Heap *heap, int node_handle)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>FMM_Heap_getHeapSize</name>
      <anchor>a11</anchor>
      <arglist>(FMM_Heap *heap)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>FMM_Heap_getHeapMemSize</name>
      <anchor>a12</anchor>
      <arglist>(FMM_Heap *heap)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FMM_Heap_printHeapData</name>
      <anchor>a13</anchor>
      <arglist>(FMM_Heap *heap)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>FMM_Macros.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/fast_marching_method/</path>
    <filename>FMM__Macros_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM_FMM_TRUE</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM_FMM_FALSE</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM_FMM_DEFAULT_UPDATE_VALUE</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM_FMM_ERR_SUCCESS</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM_FMM_ERR_FMM_DATA_CREATION_ERROR</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM_FMM_ERR_INVALID_SPATIAL_DISCRETIZATION_ORDER</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM_FMM_ABS</name>
      <anchor>a6</anchor>
      <arglist>(x)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM_FMM_IDX</name>
      <anchor>a7</anchor>
      <arglist>(idx, grid_idx, grid_dims)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM_FMM_IDX_OUT_OF_BOUNDS</name>
      <anchor>a8</anchor>
      <arglist>(result, grid_idx, grid_dims)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>LevelSetFunctionIntegrator.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/parallel/</path>
    <filename>LevelSetFunctionIntegrator_8h</filename>
    <includes id="BoundaryConditionModule_8h" name="BoundaryConditionModule.h" local="yes">BoundaryConditionModule.h</includes>
    <includes id="LevelSetFunctionIntegratorStrategy_8h" name="LevelSetFunctionIntegratorStrategy.h" local="yes">LevelSetFunctionIntegratorStrategy.h</includes>
    <includes id="OrthogonalizationAlgorithm_8h" name="OrthogonalizationAlgorithm.h" local="yes">OrthogonalizationAlgorithm.h</includes>
    <includes id="ReinitializationAlgorithm_8h" name="ReinitializationAlgorithm.h" local="yes">ReinitializationAlgorithm.h</includes>
    <includes id="LevelSetMethodToolbox_8h" name="LevelSetMethodToolbox.h" local="yes">LevelSetMethodToolbox.h</includes>
    <namespace>geom</namespace>
    <namespace>LSMLIB</namespace>
    <namespace>std</namespace>
  </compound>
  <compound kind="file">
    <name>LevelSetFunctionIntegratorStrategy.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/parallel/</path>
    <filename>LevelSetFunctionIntegratorStrategy_8h</filename>
    <includes id="LevelSetMethodToolbox_8h" name="LevelSetMethodToolbox.h" local="yes">LevelSetMethodToolbox.h</includes>
    <includes id="BoundaryConditionModule_8h" name="BoundaryConditionModule.h" local="yes">BoundaryConditionModule.h</includes>
    <namespace>LSMLIB</namespace>
    <namespace>mesh</namespace>
  </compound>
  <compound kind="file">
    <name>LevelSetMethodAlgorithm.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/parallel/</path>
    <filename>LevelSetMethodAlgorithm_8h</filename>
    <includes id="LevelSetMethodGriddingStrategy_8h" name="LevelSetMethodGriddingStrategy.h" local="yes">LevelSetMethodGriddingStrategy.h</includes>
    <includes id="LevelSetFunctionIntegratorStrategy_8h" name="LevelSetFunctionIntegratorStrategy.h" local="yes">LevelSetFunctionIntegratorStrategy.h</includes>
    <includes id="LevelSetMethodPatchStrategy_8h" name="LevelSetMethodPatchStrategy.h" local="yes">LevelSetMethodPatchStrategy.h</includes>
    <includes id="LevelSetMethodToolbox_8h" name="LevelSetMethodToolbox.h" local="yes">LevelSetMethodToolbox.h</includes>
    <includes id="LevelSetMethodVelocityFieldStrategy_8h" name="LevelSetMethodVelocityFieldStrategy.h" local="yes">LevelSetMethodVelocityFieldStrategy.h</includes>
    <includes id="FieldExtensionAlgorithm_8h" name="FieldExtensionAlgorithm.h" local="yes">FieldExtensionAlgorithm.h</includes>
    <includes id="ReinitializationAlgorithm_8h" name="ReinitializationAlgorithm.h" local="yes">ReinitializationAlgorithm.h</includes>
    <includes id="OrthogonalizationAlgorithm_8h" name="OrthogonalizationAlgorithm.h" local="yes">OrthogonalizationAlgorithm.h</includes>
    <includes id="BoundaryConditionModule_8h" name="BoundaryConditionModule.h" local="yes">BoundaryConditionModule.h</includes>
    <namespace>LSMLIB</namespace>
  </compound>
  <compound kind="file">
    <name>LevelSetMethodGriddingAlgorithm.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/parallel/</path>
    <filename>LevelSetMethodGriddingAlgorithm_8h</filename>
    <includes id="LevelSetMethodGriddingStrategy_8h" name="LevelSetMethodGriddingStrategy.h" local="yes">LevelSetMethodGriddingStrategy.h</includes>
    <includes id="LevelSetFunctionIntegratorStrategy_8h" name="LevelSetFunctionIntegratorStrategy.h" local="yes">LevelSetFunctionIntegratorStrategy.h</includes>
    <includes id="LevelSetMethodVelocityFieldStrategy_8h" name="LevelSetMethodVelocityFieldStrategy.h" local="yes">LevelSetMethodVelocityFieldStrategy.h</includes>
    <namespace>LSMLIB</namespace>
  </compound>
  <compound kind="file">
    <name>LevelSetMethodGriddingStrategy.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/parallel/</path>
    <filename>LevelSetMethodGriddingStrategy_8h</filename>
    <includes id="LevelSetMethodVelocityFieldStrategy_8h" name="LevelSetMethodVelocityFieldStrategy.h" local="yes">LevelSetMethodVelocityFieldStrategy.h</includes>
    <namespace>LSMLIB</namespace>
  </compound>
  <compound kind="file">
    <name>LevelSetMethodPatchStrategy.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/parallel/</path>
    <filename>LevelSetMethodPatchStrategy_8h</filename>
    <includes id="LevelSetFunctionIntegrator_8h" name="LevelSetFunctionIntegrator.h" local="yes">LevelSetFunctionIntegrator.h</includes>
    <includes id="LevelSetMethodVelocityFieldStrategy_8h" name="LevelSetMethodVelocityFieldStrategy.h" local="yes">LevelSetMethodVelocityFieldStrategy.h</includes>
    <namespace>LSMLIB</namespace>
  </compound>
  <compound kind="file">
    <name>LevelSetMethodToolbox.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/parallel/</path>
    <filename>LevelSetMethodToolbox_8h</filename>
    <namespace>LSMLIB</namespace>
  </compound>
  <compound kind="file">
    <name>LevelSetMethodVelocityFieldStrategy.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/parallel/</path>
    <filename>LevelSetMethodVelocityFieldStrategy_8h</filename>
    <includes id="LevelSetMethodToolbox_8h" name="LevelSetMethodToolbox.h" local="yes">LevelSetMethodToolbox.h</includes>
    <namespace>LSMLIB</namespace>
  </compound>
  <compound kind="file">
    <name>lsm_boundary_conditions.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/serial/</path>
    <filename>lsm__boundary__conditions_8h</filename>
    <includes id="lsm__grid_8h" name="lsm_grid.h" local="yes">lsm_grid.h</includes>
    <member kind="enumeration">
      <name>BOUNDARY_LOCATION_IDX</name>
      <anchor>a14</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>X_LO</name>
      <anchor>a14a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>X_HI</name>
      <anchor>a14a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>Y_LO</name>
      <anchor>a14a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>Y_HI</name>
      <anchor>a14a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>Z_LO</name>
      <anchor>a14a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>Z_HI</name>
      <anchor>a14a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>X_LO_AND_X_HI</name>
      <anchor>a14a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>Y_LO_AND_Y_HI</name>
      <anchor>a14a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>Z_LO_AND_Z_HI</name>
      <anchor>a14a8</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>ALL_BOUNDARIES</name>
      <anchor>a14a9</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>linearExtrapolationBC</name>
      <anchor>a10</anchor>
      <arglist>(LSMLIB_REAL *phi, Grid *grid, int bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>signedLinearExtrapolationBC</name>
      <anchor>a11</anchor>
      <arglist>(LSMLIB_REAL *phi, Grid *grid, int bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>copyExtrapolationBC</name>
      <anchor>a12</anchor>
      <arglist>(LSMLIB_REAL *phi, Grid *grid, int bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>homogeneousNeumannBC</name>
      <anchor>a13</anchor>
      <arglist>(LSMLIB_REAL *phi, Grid *grid, int bdry_location_idx)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_boundary_conditions1d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/boundary_conditions/</path>
    <filename>lsm__boundary__conditions1d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_LINEAR_EXTRAPOLATION</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_SIGNED_LINEAR_EXTRAPOLATION</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_COPY_EXTRAPOLATION</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_HOMOGENEOUS_NEUMANN_ENO1</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_HOMOGENEOUS_NEUMANN_ENO2</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_HOMOGENEOUS_NEUMANN_ENO3</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_HOMOGENEOUS_NEUMANN_WENO5</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_LINEAR_EXTRAPOLATION</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_SIGNED_LINEAR_EXTRAPOLATION</name>
      <anchor>a8</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_COPY_EXTRAPOLATION</name>
      <anchor>a9</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_HOMOGENEOUS_NEUMANN_ENO1</name>
      <anchor>a10</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_HOMOGENEOUS_NEUMANN_ENO2</name>
      <anchor>a11</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_HOMOGENEOUS_NEUMANN_ENO3</name>
      <anchor>a12</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_HOMOGENEOUS_NEUMANN_WENO5</name>
      <anchor>a13</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *bdry_location_idx)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_boundary_conditions2d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/boundary_conditions/</path>
    <filename>lsm__boundary__conditions2d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_LINEAR_EXTRAPOLATION</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_SIGNED_LINEAR_EXTRAPOLATION</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COPY_EXTRAPOLATION</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_HOMOGENEOUS_NEUMANN_ENO1</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_HOMOGENEOUS_NEUMANN_ENO2</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_HOMOGENEOUS_NEUMANN_ENO3</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_HOMOGENEOUS_NEUMANN_WENO5</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_LINEAR_EXTRAPOLATION</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_SIGNED_LINEAR_EXTRAPOLATION</name>
      <anchor>a8</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COPY_EXTRAPOLATION</name>
      <anchor>a9</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_HOMOGENEOUS_NEUMANN_ENO1</name>
      <anchor>a10</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_HOMOGENEOUS_NEUMANN_ENO2</name>
      <anchor>a11</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_HOMOGENEOUS_NEUMANN_ENO3</name>
      <anchor>a12</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_HOMOGENEOUS_NEUMANN_WENO5</name>
      <anchor>a13</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *bdry_location_idx)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_boundary_conditions3d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/boundary_conditions/</path>
    <filename>lsm__boundary__conditions3d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_LINEAR_EXTRAPOLATION</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_SIGNED_LINEAR_EXTRAPOLATION</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COPY_EXTRAPOLATION</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_HOMOGENEOUS_NEUMANN_ENO1</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_HOMOGENEOUS_NEUMANN_ENO2</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_HOMOGENEOUS_NEUMANN_ENO3</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_HOMOGENEOUS_NEUMANN_WENO5</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_LINEAR_EXTRAPOLATION</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *klo_gb, const  int *khi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  int *bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_SIGNED_LINEAR_EXTRAPOLATION</name>
      <anchor>a8</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *klo_gb, const  int *khi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  int *bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COPY_EXTRAPOLATION</name>
      <anchor>a9</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *klo_gb, const  int *khi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  int *bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_HOMOGENEOUS_NEUMANN_ENO1</name>
      <anchor>a10</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *klo_gb, const  int *khi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  int *bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_HOMOGENEOUS_NEUMANN_ENO2</name>
      <anchor>a11</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *klo_gb, const  int *khi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  int *bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_HOMOGENEOUS_NEUMANN_ENO3</name>
      <anchor>a12</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *klo_gb, const  int *khi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  int *bdry_location_idx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_HOMOGENEOUS_NEUMANN_WENO5</name>
      <anchor>a13</anchor>
      <arglist>(LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *klo_gb, const  int *khi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  int *bdry_location_idx)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_calculus_toolbox.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/utilities/</path>
    <filename>lsm__calculus__toolbox_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM_HEAVISIDE</name>
      <anchor>a0</anchor>
      <arglist>(x, eps)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM_DELTA_FUNCTION</name>
      <anchor>a1</anchor>
      <arglist>(x, eps)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM_HEAVISIDE_HAT</name>
      <anchor>a2</anchor>
      <arglist>(x, eps)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM_DELTA_FUNCTION_HAT</name>
      <anchor>a3</anchor>
      <arglist>(x, eps)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_calculus_toolbox2d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/utilities/</path>
    <filename>lsm__calculus__toolbox2d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_DELTA_FUNCTION_ORDER1</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_DELTA_FUNCTION_ORDER2</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_DELTA_FUNCTION_ORDER1</name>
      <anchor>a2</anchor>
      <arglist>(const  LSMLIB_REAL *phi, const  LSMLIB_REAL *delta, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  LSMLIB_REAL *norm_phi_x, const  LSMLIB_REAL *norm_phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_DELTA_FUNCTION_ORDER2</name>
      <anchor>a3</anchor>
      <arglist>(const  LSMLIB_REAL *phi, const  LSMLIB_REAL *delta, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  LSMLIB_REAL *norm_phi_x, const  LSMLIB_REAL *norm_phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_calculus_toolbox2d_local.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/utilities/</path>
    <filename>lsm__calculus__toolbox2d__local_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_DELTA_FUNCTION_ORDER1_LOCAL</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_DELTA_FUNCTION_ORDER2_LOCAL</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_DELTA_FUNCTION_ORDER1_LOCAL</name>
      <anchor>a2</anchor>
      <arglist>(const  LSMLIB_REAL *phi, const  LSMLIB_REAL *delta, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  LSMLIB_REAL *norm_phi_x, const  LSMLIB_REAL *norm_phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_DELTA_FUNCTION_ORDER2_LOCAL</name>
      <anchor>a3</anchor>
      <arglist>(const  LSMLIB_REAL *phi, const  LSMLIB_REAL *delta, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  LSMLIB_REAL *norm_phi_x, const  LSMLIB_REAL *norm_phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_calculus_toolbox3d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/utilities/</path>
    <filename>lsm__calculus__toolbox3d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_DELTA_FUNCTION_ORDER1</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_DELTA_FUNCTION_ORDER1</name>
      <anchor>a1</anchor>
      <arglist>(const  LSMLIB_REAL *phi, const  LSMLIB_REAL *delta, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *klo_gb, const  int *khi_gb, const  LSMLIB_REAL *norm_phi_x, const  LSMLIB_REAL *norm_phi_y, const  LSMLIB_REAL *norm_phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_curvature2d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/geometry/</path>
    <filename>lsm__curvature2d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_MEAN_CURVATURE_ORDER2</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_MEAN_CURVATURE_9STENCIL</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_MEAN_CURVATURE_SGN_DIST</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_MEAN_CURVATURE_ORDER4</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_MEAN_CURVATURE_ORDER2</name>
      <anchor>a4</anchor>
      <arglist>(LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_MEAN_CURVATURE_9STENCIL</name>
      <anchor>a5</anchor>
      <arglist>(LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_MEAN_CURVATURE_SGN_DIST</name>
      <anchor>a6</anchor>
      <arglist>(LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_MEAN_CURVATURE_ORDER4</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_curvature2d_local.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/geometry/</path>
    <filename>lsm__curvature2d__local_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_MEAN_CURVATURE_ORDER2_LOCAL</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_MEAN_CURVATURE_9STENCIL_LOCAL</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_MEAN_CURVATURE_SGN_DIST_LOCAL</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_MEAN_CURVATURE_ORDER4_LOCAL</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_MEAN_CURVATURE_ORDER2_LOCAL</name>
      <anchor>a4</anchor>
      <arglist>(LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_MEAN_CURVATURE_9STENCIL_LOCAL</name>
      <anchor>a5</anchor>
      <arglist>(LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_MEAN_CURVATURE_SGN_DIST_LOCAL</name>
      <anchor>a6</anchor>
      <arglist>(LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_MEAN_CURVATURE_ORDER4_LOCAL</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_curvature3d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/geometry/</path>
    <filename>lsm__curvature3d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_MEAN_CURVATURE_ORDER2</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_GAUSSIAN_CURVATURE_ORDER2</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_MEAN_CURVATURE_ORDER4</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_GAUSSIAN_CURVATURE_ORDER4</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_MEAN_CURVATURE_ORDER2</name>
      <anchor>a4</anchor>
      <arglist>(LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  int *klo_kappa_gb, const  int *khi_kappa_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  int *ilo_kappa_fb, const  int *ihi_kappa_fb, const  int *jlo_kappa_fb, const  int *jhi_kappa_fb, const  int *klo_kappa_fb, const  int *khi_kappa_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_GAUSSIAN_CURVATURE_ORDER2</name>
      <anchor>a5</anchor>
      <arglist>(LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  int *klo_kappa_gb, const  int *khi_kappa_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  int *ilo_kappa_fb, const  int *ihi_kappa_fb, const  int *jlo_kappa_fb, const  int *jhi_kappa_fb, const  int *klo_kappa_fb, const  int *khi_kappa_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_MEAN_CURVATURE_ORDER4</name>
      <anchor>a6</anchor>
      <arglist>(LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  int *klo_kappa_gb, const  int *khi_kappa_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  int *ilo_kappa_fb, const  int *ihi_kappa_fb, const  int *jlo_kappa_fb, const  int *jhi_kappa_fb, const  int *klo_kappa_fb, const  int *khi_kappa_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_GAUSSIAN_CURVATURE_ORDER4</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  int *klo_kappa_gb, const  int *khi_kappa_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  int *ilo_kappa_fb, const  int *ihi_kappa_fb, const  int *jlo_kappa_fb, const  int *jhi_kappa_fb, const  int *klo_kappa_fb, const  int *khi_kappa_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_curvature3d_local.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/geometry/</path>
    <filename>lsm__curvature3d__local_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_MEAN_CURVATURE_ORDER2_LOCAL</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_GAUSSIAN_CURVATURE_ORDER2_LOCAL</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_MEAN_CURVATURE_ORDER4_LOCAL</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_GAUSSIAN_CURVATURE_ORDER4_LOCAL</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_MEAN_CURVATURE_ORDER2_LOCAL</name>
      <anchor>a4</anchor>
      <arglist>(LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  int *klo_kappa_gb, const  int *khi_kappa_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_GAUSSIAN_CURVATURE_ORDER2_LOCAL</name>
      <anchor>a5</anchor>
      <arglist>(LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  int *klo_kappa_gb, const  int *khi_kappa_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_MEAN_CURVATURE_ORDER4_LOCAL</name>
      <anchor>a6</anchor>
      <arglist>(LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  int *klo_kappa_gb, const  int *khi_kappa_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_GAUSSIAN_CURVATURE_ORDER4_LOCAL</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  int *klo_kappa_gb, const  int *khi_kappa_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_data_arrays.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/serial/</path>
    <filename>lsm__data__arrays_8h</filename>
    <includes id="lsm__grid_8h" name="lsm_grid.h" local="yes">lsm_grid.h</includes>
    <includes id="lsm__file_8h" name="lsm_file.h" local="yes">lsm_file.h</includes>
    <class kind="struct">_LSM_DataArrays</class>
    <member kind="typedef">
      <type>_LSM_DataArrays</type>
      <name>LSM_DataArrays</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>LSM_DataArrays *</type>
      <name>allocateLSMDataArrays</name>
      <anchor>a1</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>allocateMemoryForLSMDataArrays</name>
      <anchor>a2</anchor>
      <arglist>(LSM_DataArrays *lsm_data_arrays, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>freeMemoryForLSMDataArrays</name>
      <anchor>a3</anchor>
      <arglist>(LSM_DataArrays *lsm_arrays)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>destroyLSMDataArrays</name>
      <anchor>a4</anchor>
      <arglist>(LSM_DataArrays *lsm_data_arrays)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>writeDataArray</name>
      <anchor>a5</anchor>
      <arglist>(LSMLIB_REAL *data, Grid *grid, char *file_name, int zip_status)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>writeDataArrayNoGrid</name>
      <anchor>a6</anchor>
      <arglist>(LSMLIB_REAL *data, int *n, char *file_name, int zip_status)</arglist>
    </member>
    <member kind="function">
      <type>LSMLIB_REAL *</type>
      <name>readDataArray</name>
      <anchor>a7</anchor>
      <arglist>(int *grid_dims, char *file_name)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>writeDataArray1d</name>
      <anchor>a8</anchor>
      <arglist>(LSMLIB_REAL *data, int num_gridpts, char *file_name, int zip_status)</arglist>
    </member>
    <member kind="function">
      <type>LSMLIB_REAL *</type>
      <name>readDataArray1d</name>
      <anchor>a9</anchor>
      <arglist>(int *num_elements, char *file_name)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_fast_marching_method.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/serial/</path>
    <filename>lsm__fast__marching__method_8h</filename>
    <includes id="FMM__Core_8h" name="FMM_Core.h" local="yes">FMM_Core.h</includes>
    <member kind="function">
      <type>int</type>
      <name>computeExtensionFields2d</name>
      <anchor>a0</anchor>
      <arglist>(LSMLIB_REAL *distance_function, LSMLIB_REAL **extension_fields, LSMLIB_REAL *phi, LSMLIB_REAL *mark, LSMLIB_REAL **source_fields, int num_extension_fields, int spatial_discretization_order, int *grid_dims, LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>computeDistanceFunction2d</name>
      <anchor>a1</anchor>
      <arglist>(LSMLIB_REAL *distance_function, LSMLIB_REAL *phi, LSMLIB_REAL *mark, int spatial_discretization_order, int *grid_dims, LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>solveEikonalEquation2d</name>
      <anchor>a2</anchor>
      <arglist>(LSMLIB_REAL *phi, LSMLIB_REAL *speed, LSMLIB_REAL *mask, int spatial_discretization_order, int *grid_dims, LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>computeExtensionFields3d</name>
      <anchor>a3</anchor>
      <arglist>(LSMLIB_REAL *distance_function, LSMLIB_REAL **extension_fields, LSMLIB_REAL *phi, LSMLIB_REAL *mask, LSMLIB_REAL **source_fields, int num_extension_fields, int spatial_discretization_order, int *grid_dims, LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>computeDistanceFunction3d</name>
      <anchor>a4</anchor>
      <arglist>(LSMLIB_REAL *distance_function, LSMLIB_REAL *phi, LSMLIB_REAL *mask, int spatial_discretization_order, int *grid_dims, LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>solveEikonalEquation3d</name>
      <anchor>a5</anchor>
      <arglist>(LSMLIB_REAL *phi, LSMLIB_REAL *speed, LSMLIB_REAL *mask, int spatial_discretization_order, int *grid_dims, LSMLIB_REAL *dx)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_field_extension1d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/field_extension/</path>
    <filename>lsm__field__extension1d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_COMPUTE_FIELD_EXTENSION_EQN_RHS</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_COMPUTE_FIELD_EXTENSION_EQN_RHS</name>
      <anchor>a1</anchor>
      <arglist>(LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  LSMLIB_REAL *S, const  int *ilo_S_gb, const  int *ihi_S_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  LSMLIB_REAL *S_x_upwind, const  int *ilo_grad_S_upwind_gb, const  int *ihi_grad_S_upwind_gb, const  LSMLIB_REAL *signed_normal_x, const  int *ilo_signed_normal_gb, const  int *ihi_signed_normal_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dx)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_field_extension2d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/field_extension/</path>
    <filename>lsm__field__extension2d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_FIELD_EXTENSION_EQN_RHS</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_FIELD_EXTENSION_EQN_RHS</name>
      <anchor>a1</anchor>
      <arglist>(LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  LSMLIB_REAL *S, const  int *ilo_S_gb, const  int *ihi_S_gb, const  int *jlo_S_gb, const  int *jhi_S_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *S_x_upwind, const  LSMLIB_REAL *S_y_upwind, const  int *ilo_grad_S_upwind_gb, const  int *ihi_grad_S_upwind_gb, const  int *jlo_grad_S_upwind_gb, const  int *jhi_grad_S_upwind_gb, const  LSMLIB_REAL *signed_normal_x, const  LSMLIB_REAL *signed_normal_y, const  int *ilo_signed_normal_gb, const  int *ihi_signed_normal_gb, const  int *jlo_signed_normal_gb, const  int *jhi_signed_normal_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_field_extension3d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/field_extension/</path>
    <filename>lsm__field__extension3d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_FIELD_EXTENSION_EQN_RHS</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_FIELD_EXTENSION_EQN_RHS</name>
      <anchor>a1</anchor>
      <arglist>(LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  LSMLIB_REAL *S, const  int *ilo_S_gb, const  int *ihi_S_gb, const  int *jlo_S_gb, const  int *jhi_S_gb, const  int *klo_S_gb, const  int *khi_S_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *S_x_upwind, const  LSMLIB_REAL *S_y_upwind, const  LSMLIB_REAL *S_z_upwind, const  int *ilo_grad_S_upwind_gb, const  int *ihi_grad_S_upwind_gb, const  int *jlo_grad_S_upwind_gb, const  int *jhi_grad_S_upwind_gb, const  int *klo_grad_S_upwind_gb, const  int *khi_grad_S_upwind_gb, const  LSMLIB_REAL *signed_normal_x, const  LSMLIB_REAL *signed_normal_y, const  LSMLIB_REAL *signed_normal_z, const  int *ilo_signed_normal_gb, const  int *ihi_signed_normal_gb, const  int *jlo_signed_normal_gb, const  int *jhi_signed_normal_gb, const  int *klo_signed_normal_gb, const  int *khi_signed_normal_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_file.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/serial/</path>
    <filename>lsm__file_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>NO_ZIP</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>GZIP</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>BZIP2</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>checkUnzipFile</name>
      <anchor>a3</anchor>
      <arglist>(char *file_name, int *pzip_status, char **pfile_base)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>zipFile</name>
      <anchor>a4</anchor>
      <arglist>(char *file_base, int zip_status)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_geometry1d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/geometry/</path>
    <filename>lsm__geometry1d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_COMPUTE_UNIT_NORMAL</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_COMPUTE_SIGNED_UNIT_NORMAL</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_LENGTH_REGION_PHI_LESS_THAN_ZERO</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_LENGTH_REGION_PHI_GREATER_THAN_ZERO</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_SIZE_ZERO_LEVEL_SET</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_LENGTH_REGION_PHI_LESS_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_LENGTH_REGION_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_SIZE_ZERO_LEVEL_SET_CONTROL_VOLUME</name>
      <anchor>a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_COMPUTE_UNIT_NORMAL</name>
      <anchor>a8</anchor>
      <arglist>(LSMLIB_REAL *normal, const  int *ilo_normal_gb, const  int *ihi_normal_gb, const  LSMLIB_REAL *phi_x, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *ilo_fb, const  int *ihi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_COMPUTE_SIGNED_UNIT_NORMAL</name>
      <anchor>a9</anchor>
      <arglist>(LSMLIB_REAL *normal, const  int *ilo_normal_gb, const  int *ihi_normal_gb, const  LSMLIB_REAL *phi_x, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_LENGTH_REGION_PHI_LESS_THAN_ZERO</name>
      <anchor>a10</anchor>
      <arglist>(LSMLIB_REAL *length, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_LENGTH_REGION_PHI_GREATER_THAN_ZERO</name>
      <anchor>a11</anchor>
      <arglist>(LSMLIB_REAL *volume, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_SIZE_ZERO_LEVEL_SET</name>
      <anchor>a12</anchor>
      <arglist>(LSMLIB_REAL *size, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  LSMLIB_REAL *phi_x, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_LENGTH_REGION_PHI_LESS_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a13</anchor>
      <arglist>(LSMLIB_REAL *length, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_LENGTH_REGION_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a14</anchor>
      <arglist>(LSMLIB_REAL *volume, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_SIZE_ZERO_LEVEL_SET_CONTROL_VOLUME</name>
      <anchor>a15</anchor>
      <arglist>(LSMLIB_REAL *size, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  LSMLIB_REAL *phi_x, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_geometry2d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/geometry/</path>
    <filename>lsm__geometry2d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_UNIT_NORMAL</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_SIGNED_UNIT_NORMAL</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_AREA_REGION_PHI_LESS_THAN_ZERO</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_AREA_REGION_PHI_GREATER_THAN_ZERO</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_PERIMETER_ZERO_LEVEL_SET</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_PERIMETER_ZERO_LEVEL_SET_DELTA</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_AREA_REGION_PHI_LESS_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_AREA_REGION_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_PERIMETER_ZERO_LEVEL_SET_CONTROL_VOLUME</name>
      <anchor>a8</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_PERIMETER_ZERO_LEVEL_SET_DELTA_CONTROL_VOLUME</name>
      <anchor>a9</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_UNIT_NORMAL</name>
      <anchor>a10</anchor>
      <arglist>(LSMLIB_REAL *normal_x, LSMLIB_REAL *normal_y, const  int *ilo_normal_gb, const  int *ihi_normal_gb, const  int *jlo_normal_gb, const  int *jhi_normal_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_SIGNED_UNIT_NORMAL</name>
      <anchor>a11</anchor>
      <arglist>(LSMLIB_REAL *normal_x, LSMLIB_REAL *normal_y, const  int *ilo_normal_gb, const  int *ihi_normal_gb, const  int *jlo_normal_gb, const  int *jhi_normal_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_AREA_REGION_PHI_LESS_THAN_ZERO</name>
      <anchor>a12</anchor>
      <arglist>(LSMLIB_REAL *area, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_AREA_REGION_PHI_GREATER_THAN_ZERO</name>
      <anchor>a13</anchor>
      <arglist>(LSMLIB_REAL *area, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_PERIMETER_ZERO_LEVEL_SET</name>
      <anchor>a14</anchor>
      <arglist>(LSMLIB_REAL *perimeter, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_PERIMETER_ZERO_LEVEL_SET_DELTA</name>
      <anchor>a15</anchor>
      <arglist>(LSMLIB_REAL *perimeter, const  LSMLIB_REAL *delta_phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_AREA_REGION_PHI_LESS_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a16</anchor>
      <arglist>(LSMLIB_REAL *area, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_AREA_REGION_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a17</anchor>
      <arglist>(LSMLIB_REAL *area, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_PERIMETER_ZERO_LEVEL_SET_CONTROL_VOLUME</name>
      <anchor>a18</anchor>
      <arglist>(LSMLIB_REAL *perimeter, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_PERIMETER_ZERO_LEVEL_SET_DELTA_CONTROL_VOLUME</name>
      <anchor>a19</anchor>
      <arglist>(LSMLIB_REAL *perimeter, const  LSMLIB_REAL *delta_phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_geometry2d_local.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/geometry/</path>
    <filename>lsm__geometry2d__local_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_UNIT_NORMAL_LOCAL</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_SIGNED_UNIT_NORMAL_LOCAL</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_PERIMETER_ZERO_LEVEL_SET_DELTA_CONTROL_VOLUME_LOCAL</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_UNIT_NORMAL_LOCAL</name>
      <anchor>a3</anchor>
      <arglist>(LSMLIB_REAL *normal_x, LSMLIB_REAL *normal_y, const  int *ilo_normal_gb, const  int *ihi_normal_gb, const  int *jlo_normal_gb, const  int *jhi_normal_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_SIGNED_UNIT_NORMAL_LOCAL</name>
      <anchor>a4</anchor>
      <arglist>(LSMLIB_REAL *normal_x, LSMLIB_REAL *normal_y, const  int *ilo_normal_gb, const  int *ihi_normal_gb, const  int *jlo_normal_gb, const  int *jhi_normal_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_PERIMETER_ZERO_LEVEL_SET_DELTA_CONTROL_VOLUME_LOCAL</name>
      <anchor>a5</anchor>
      <arglist>(LSMLIB_REAL *perimeter, const  LSMLIB_REAL *delta_phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *control_vol_sgn, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_geometry3d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/geometry/</path>
    <filename>lsm__geometry3d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_UNIT_NORMAL</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_SIGNED_UNIT_NORMAL</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_VOLUME_REGION_PHI_LESS_THAN_ZERO</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_VOLUME_REGION_PHI_GREATER_THAN_ZERO</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_SURFACE_AREA_ZERO_LEVEL_SET</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_SURFACE_AREA_ZERO_LEVEL_SET_DELTA</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_VOLUME_REGION_PHI_LESS_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_VOLUME_REGION_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_SURFACE_AREA_ZERO_LEVEL_SET_CONTROL_VOLUME</name>
      <anchor>a8</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_SURFACE_AREA_ZERO_LEVEL_SET_DELTA_CONTROL_VOLUME</name>
      <anchor>a9</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_UNIT_NORMAL</name>
      <anchor>a10</anchor>
      <arglist>(LSMLIB_REAL *normal_x, LSMLIB_REAL *normal_y, LSMLIB_REAL *normal_z, const  int *ilo_normal_gb, const  int *ihi_normal_gb, const  int *jlo_normal_gb, const  int *jhi_normal_gb, const  int *klo_normal_gb, const  int *khi_normal_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_SIGNED_UNIT_NORMAL</name>
      <anchor>a11</anchor>
      <arglist>(LSMLIB_REAL *normal_x, LSMLIB_REAL *normal_y, LSMLIB_REAL *normal_z, const  int *ilo_normal_gb, const  int *ihi_normal_gb, const  int *jlo_normal_gb, const  int *jhi_normal_gb, const  int *klo_normal_gb, const  int *khi_normal_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_VOLUME_REGION_PHI_LESS_THAN_ZERO</name>
      <anchor>a12</anchor>
      <arglist>(LSMLIB_REAL *volume, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_VOLUME_REGION_PHI_GREATER_THAN_ZERO</name>
      <anchor>a13</anchor>
      <arglist>(LSMLIB_REAL *volume, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_SURFACE_AREA_ZERO_LEVEL_SET</name>
      <anchor>a14</anchor>
      <arglist>(LSMLIB_REAL *surface_area, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_SURFACE_AREA_ZERO_LEVEL_SET_DELTA</name>
      <anchor>a15</anchor>
      <arglist>(LSMLIB_REAL *surface_area, const  LSMLIB_REAL *delta_phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_VOLUME_REGION_PHI_LESS_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a16</anchor>
      <arglist>(LSMLIB_REAL *volume, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *klo_control_vol_gb, const  int *khi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_VOLUME_REGION_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a17</anchor>
      <arglist>(LSMLIB_REAL *volume, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *klo_control_vol_gb, const  int *khi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_SURFACE_AREA_ZERO_LEVEL_SET_CONTROL_VOLUME</name>
      <anchor>a18</anchor>
      <arglist>(LSMLIB_REAL *surface_area, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *klo_control_vol_gb, const  int *khi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_SURFACE_AREA_ZERO_LEVEL_SET_DELTA_CONTROL_VOLUME</name>
      <anchor>a19</anchor>
      <arglist>(LSMLIB_REAL *surface_area, const  LSMLIB_REAL *delta_phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *klo_control_vol_gb, const  int *khi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>LSM3D_findLineInTetrahedron</name>
      <anchor>a20</anchor>
      <arglist>(LSMLIB_REAL *endpt1, LSMLIB_REAL *endpt2, const  LSMLIB_REAL *x1, const  LSMLIB_REAL *x2, const  LSMLIB_REAL *x3, const  LSMLIB_REAL *x4, const  LSMLIB_REAL *phi, const  LSMLIB_REAL *psi)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_grid.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/serial/</path>
    <filename>lsm__grid_8h</filename>
    <class kind="struct">_Grid</class>
    <member kind="function">
      <type>Grid *</type>
      <name>createGridSetDx</name>
      <anchor>a5</anchor>
      <arglist>(int num_dims, LSMLIB_REAL dx, LSMLIB_REAL *x_lo, LSMLIB_REAL *x_hi, LSMLIB_SPATIAL_DERIVATIVE_ACCURACY_TYPE accuracy)</arglist>
    </member>
    <member kind="function">
      <type>Grid *</type>
      <name>createGridSetDxDyDz</name>
      <anchor>a6</anchor>
      <arglist>(int num_dims, LSMLIB_REAL dx, LSMLIB_REAL dy, LSMLIB_REAL dz, LSMLIB_REAL *x_lo, LSMLIB_REAL *x_hi, LSMLIB_SPATIAL_DERIVATIVE_ACCURACY_TYPE accuracy)</arglist>
    </member>
    <member kind="function">
      <type>Grid *</type>
      <name>createGridSetGridDims</name>
      <anchor>a7</anchor>
      <arglist>(int num_dims, int *grid_dims, LSMLIB_REAL *x_lo, LSMLIB_REAL *x_hi, LSMLIB_SPATIAL_DERIVATIVE_ACCURACY_TYPE accuracy)</arglist>
    </member>
    <member kind="function">
      <type>Grid *</type>
      <name>copyGrid</name>
      <anchor>a8</anchor>
      <arglist>(Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>destroyGrid</name>
      <anchor>a9</anchor>
      <arglist>(Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>printGrid</name>
      <anchor>a10</anchor>
      <arglist>(Grid *grid, FILE *fp)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>writeGridToAsciiFile</name>
      <anchor>a11</anchor>
      <arglist>(Grid *grid, char *file_name, int zip_status)</arglist>
    </member>
    <member kind="function">
      <type>Grid *</type>
      <name>readGridFromAsciiFile</name>
      <anchor>a12</anchor>
      <arglist>(char *file_name)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>writeGridToBinaryFile</name>
      <anchor>a13</anchor>
      <arglist>(Grid *grid, char *file_name, int zip_status)</arglist>
    </member>
    <member kind="function">
      <type>Grid *</type>
      <name>readGridFromBinaryFile</name>
      <anchor>a14</anchor>
      <arglist>(char *file_name)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setIndexSpaceLimits</name>
      <anchor>a15</anchor>
      <arglist>(LSMLIB_SPATIAL_DERIVATIVE_ACCURACY_TYPE accuracy, Grid *grid)</arglist>
    </member>
    <member kind="typedef">
      <type>_Grid</type>
      <name>Grid</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumeration">
      <name>LSMLIB_SPATIAL_DERIVATIVE_ACCURACY_TYPE</name>
      <anchor>a16</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>LOW</name>
      <anchor>a16a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>MEDIUM</name>
      <anchor>a16a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>HIGH</name>
      <anchor>a16a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>VERY_HIGH</name>
      <anchor>a16a4</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_initialization2d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/serial/</path>
    <filename>lsm__initialization2d_8h</filename>
    <includes id="lsm__grid_8h" name="lsm_grid.h" local="yes">lsm_grid.h</includes>
    <member kind="function">
      <type>void</type>
      <name>createLine</name>
      <anchor>a0</anchor>
      <arglist>(LSMLIB_REAL *phi, LSMLIB_REAL normal_x, LSMLIB_REAL normal_y, LSMLIB_REAL point_x, LSMLIB_REAL point_y, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createIntersectionOfHalfSpaces2d</name>
      <anchor>a1</anchor>
      <arglist>(LSMLIB_REAL *phi, int num_half_spaces, LSMLIB_REAL *normal_x, LSMLIB_REAL *normal_y, LSMLIB_REAL *point_x, LSMLIB_REAL *point_y, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createPolyhedron2d</name>
      <anchor>a2</anchor>
      <arglist>(LSMLIB_REAL *phi, int num_sides, LSMLIB_REAL *normal_x, LSMLIB_REAL *normal_y, LSMLIB_REAL *point_x, LSMLIB_REAL *point_y, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createIntersectionOfPolyhedra2d</name>
      <anchor>a3</anchor>
      <arglist>(LSMLIB_REAL *phi, int num_polyhedra, int *idx_start, int *idx_end, LSMLIB_REAL *normal_x, LSMLIB_REAL *normal_y, LSMLIB_REAL *point_x, LSMLIB_REAL *point_y, int *inside_flag, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createCircle</name>
      <anchor>a4</anchor>
      <arglist>(LSMLIB_REAL *phi, LSMLIB_REAL center_x, LSMLIB_REAL center_y, LSMLIB_REAL radius, int inside_flag, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createIntersectionOfCircles</name>
      <anchor>a5</anchor>
      <arglist>(LSMLIB_REAL *phi, int num_circles, LSMLIB_REAL *center_x, LSMLIB_REAL *center_y, LSMLIB_REAL *radius, int *inside_flag, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createRectangle</name>
      <anchor>a6</anchor>
      <arglist>(LSMLIB_REAL *phi, LSMLIB_REAL corner_x, LSMLIB_REAL corner_y, LSMLIB_REAL side_length_x, LSMLIB_REAL side_length_y, int inside_flag, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createIntersectionOfRectangles</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *phi, int num_rectangles, LSMLIB_REAL *corner_x, LSMLIB_REAL *corner_y, LSMLIB_REAL *side_length_x, LSMLIB_REAL *side_length_y, int *inside_flag, Grid *grid)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_initialization3d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/serial/</path>
    <filename>lsm__initialization3d_8h</filename>
    <includes id="lsm__grid_8h" name="lsm_grid.h" local="yes">lsm_grid.h</includes>
    <member kind="function">
      <type>void</type>
      <name>createPlane</name>
      <anchor>a0</anchor>
      <arglist>(LSMLIB_REAL *phi, LSMLIB_REAL normal_x, LSMLIB_REAL normal_y, LSMLIB_REAL normal_z, LSMLIB_REAL point_x, LSMLIB_REAL point_y, LSMLIB_REAL point_z, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createIntersectionOfHalfSpaces3d</name>
      <anchor>a1</anchor>
      <arglist>(LSMLIB_REAL *phi, int num_half_spaces, LSMLIB_REAL *normal_x, LSMLIB_REAL *normal_y, LSMLIB_REAL *normal_z, LSMLIB_REAL *point_x, LSMLIB_REAL *point_y, LSMLIB_REAL *point_z, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createPolyhedron3d</name>
      <anchor>a2</anchor>
      <arglist>(LSMLIB_REAL *phi, int num_faces, LSMLIB_REAL *normal_x, LSMLIB_REAL *normal_y, LSMLIB_REAL *normal_z, LSMLIB_REAL *point_x, LSMLIB_REAL *point_y, LSMLIB_REAL *point_z, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createSphere</name>
      <anchor>a3</anchor>
      <arglist>(LSMLIB_REAL *phi, LSMLIB_REAL center_x, LSMLIB_REAL center_y, LSMLIB_REAL center_z, LSMLIB_REAL radius, int inside_flag, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createIntersectionOfSpheres</name>
      <anchor>a4</anchor>
      <arglist>(LSMLIB_REAL *phi, int num_spheres, LSMLIB_REAL *center_x, LSMLIB_REAL *center_y, LSMLIB_REAL *center_z, LSMLIB_REAL *radius, int *inside_flag, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createCylinder</name>
      <anchor>a5</anchor>
      <arglist>(LSMLIB_REAL *phi, LSMLIB_REAL tangent_x, LSMLIB_REAL tangent_y, LSMLIB_REAL tangent_z, LSMLIB_REAL point_x, LSMLIB_REAL point_y, LSMLIB_REAL point_z, LSMLIB_REAL radius, int inside_flag, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createIntersectionOfCylinders</name>
      <anchor>a6</anchor>
      <arglist>(LSMLIB_REAL *phi, int num_cylinders, LSMLIB_REAL *tangent_x, LSMLIB_REAL *tangent_y, LSMLIB_REAL *tangent_z, LSMLIB_REAL *point_x, LSMLIB_REAL *point_y, LSMLIB_REAL *point_z, LSMLIB_REAL *radius, int *inside_flag, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createHyperboloid</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *phi, LSMLIB_REAL tangent_x, LSMLIB_REAL tangent_y, LSMLIB_REAL tangent_z, LSMLIB_REAL center_x, LSMLIB_REAL center_y, LSMLIB_REAL center_z, LSMLIB_REAL alpha, LSMLIB_REAL beta, int inside_flag, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createIntersectionOfHyperboloids</name>
      <anchor>a8</anchor>
      <arglist>(LSMLIB_REAL *phi, int num_hyperboloids, LSMLIB_REAL *tangent_x, LSMLIB_REAL *tangent_y, LSMLIB_REAL *tangent_z, LSMLIB_REAL *center_x, LSMLIB_REAL *center_y, LSMLIB_REAL *center_z, LSMLIB_REAL *alpha, LSMLIB_REAL *beta, int *inside_flag, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createCone</name>
      <anchor>a9</anchor>
      <arglist>(LSMLIB_REAL *phi, LSMLIB_REAL tangent_x, LSMLIB_REAL tangent_y, LSMLIB_REAL tangent_z, LSMLIB_REAL center_x, LSMLIB_REAL center_y, LSMLIB_REAL center_z, LSMLIB_REAL alpha, LSMLIB_REAL beta, int inside_flag, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createIntersectionOfCones</name>
      <anchor>a10</anchor>
      <arglist>(LSMLIB_REAL *phi, int num_cones, LSMLIB_REAL *tangent_x, LSMLIB_REAL *tangent_y, LSMLIB_REAL *tangent_z, LSMLIB_REAL *center_x, LSMLIB_REAL *center_y, LSMLIB_REAL *center_z, LSMLIB_REAL *alpha, LSMLIB_REAL *beta, int *inside_flag, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createBox</name>
      <anchor>a11</anchor>
      <arglist>(LSMLIB_REAL *phi, LSMLIB_REAL corner_x, LSMLIB_REAL corner_y, LSMLIB_REAL corner_z, LSMLIB_REAL side_length_x, LSMLIB_REAL side_length_y, LSMLIB_REAL side_length_z, int inside_flag, Grid *grid)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>createIntersectionOfBoxes</name>
      <anchor>a12</anchor>
      <arglist>(LSMLIB_REAL *phi, int num_cuboids, LSMLIB_REAL *corner_x, LSMLIB_REAL *corner_y, LSMLIB_REAL *corner_z, LSMLIB_REAL *side_length_x, LSMLIB_REAL *side_length_y, LSMLIB_REAL *side_length_z, int *inside_flag, Grid *grid)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_level_set_evolution1d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/level_set_evolution/</path>
    <filename>lsm__level__set__evolution1d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_ZERO_OUT_LEVEL_SET_EQN_RHS</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_ADD_ADVECTION_TERM_TO_LSE_RHS</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_ZERO_OUT_LEVEL_SET_EQN_RHS</name>
      <anchor>a3</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_ADD_ADVECTION_TERM_TO_LSE_RHS</name>
      <anchor>a4</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  LSMLIB_REAL *phi_x, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  LSMLIB_REAL *vel_x, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *ilo_fb, const  int *ihi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS</name>
      <anchor>a5</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  LSMLIB_REAL *phi_x_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  LSMLIB_REAL *vel_n, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *ilo_fb, const  int *ihi_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_level_set_evolution2d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/level_set_evolution/</path>
    <filename>lsm__level__set__evolution2d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_ZERO_OUT_LEVEL_SET_EQN_RHS</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_ADD_ADVECTION_TERM_TO_LSE_RHS</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_ADD_CONST_CURV_TERM_TO_LSE_RHS</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_ADD_CONST_PRECOMPUTED_CURV_TERM_TO_LSE_RHS</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_ADD_EXTERNAL_AND_NORMAL_VEL_TERM_TO_LSE_RHS</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_ZERO_OUT_LEVEL_SET_EQN_RHS</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_ADD_ADVECTION_TERM_TO_LSE_RHS</name>
      <anchor>a8</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS</name>
      <anchor>a9</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *vel_n, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS</name>
      <anchor>a10</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *vel_n, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_ADD_CONST_CURV_TERM_TO_LSE_RHS</name>
      <anchor>a11</anchor>
      <arglist>(const  LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *phi_xx, const  LSMLIB_REAL *phi_xy, const  LSMLIB_REAL *phi_yy, const  int *ilo_grad2_phi_gb, const  int *ihi_grad2_phi_gb, const  int *jlo_grad2_phi_gb, const  int *jhi_grad2_phi_gb, const  LSMLIB_REAL *b, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_ADD_CONST_PRECOMPUTED_CURV_TERM_TO_LSE_RHS</name>
      <anchor>a12</anchor>
      <arglist>(const  LSMLIB_REAL *lse_rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  LSMLIB_REAL *grad_mag_phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *b, const  int *ilo_rhs_fb, const  int *ihi_rhs_fb, const  int *jlo_rhs_fb, const  int *jhi_rhs_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_ADD_EXTERNAL_AND_NORMAL_VEL_TERM_TO_LSE_RHS</name>
      <anchor>a13</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *vel_n, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *ilo_rhs_fb, const  int *ihi_rhs_fb, const  int *jlo_rhs_fb, const  int *jhi_rhs_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_level_set_evolution2d_local.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/level_set_evolution/</path>
    <filename>lsm__level__set__evolution2d__local_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_ZERO_OUT_LEVEL_SET_EQN_RHS_LOCAL</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_ADD_ADVECTION_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_ADD_CONST_CURV_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_ADD_CONST_PRECOMPUTED_CURV_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_ADD_EXTERNAL_AND_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_ZERO_OUT_LEVEL_SET_EQN_RHS_LOCAL</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_ADD_ADVECTION_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a8</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a9</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *vel_n, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a10</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *vel_n, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_ADD_CONST_CURV_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a11</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *phi_xx, const  LSMLIB_REAL *phi_xy, const  LSMLIB_REAL *phi_yy, const  int *ilo_grad2_phi_gb, const  int *ihi_grad2_phi_gb, const  int *jlo_grad2_phi_gb, const  int *jhi_grad2_phi_gb, const  LSMLIB_REAL *b, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_ADD_CONST_PRECOMPUTED_CURV_TERM_TO_LSERHS_LOCAL</name>
      <anchor>a12</anchor>
      <arglist>(const  LSMLIB_REAL *lse_rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  LSMLIB_REAL *grad_mag_phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *b, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_ADD_EXTERNAL_AND_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a13</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *vel_n, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_level_set_evolution3d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/level_set_evolution/</path>
    <filename>lsm__level__set__evolution3d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_ZERO_OUT_LEVEL_SET_EQN_RHS</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_ADD_ADVECTION_TERM_TO_LSE_RHS</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_ADD_CONST_CURV_TERM_TO_LSE_RHS</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_ADD_CONST_PRECOMPUTED_CURV_TERM_TO_LSE_RHS</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_ADD_EXTERNAL_AND_NORMAL_VEL_TERM_TO_LSE_RHS</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_ZERO_OUT_LEVEL_SET_EQN_RHS</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb, const  int *klo_lse_rhs_gb, const  int *khi_lse_rhs_gb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_ADD_ADVECTION_TERM_TO_LSE_RHS</name>
      <anchor>a8</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb, const  int *klo_lse_rhs_gb, const  int *khi_lse_rhs_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  LSMLIB_REAL *vel_z, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *klo_vel_gb, const  int *khi_vel_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS</name>
      <anchor>a9</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb, const  int *klo_lse_rhs_gb, const  int *khi_lse_rhs_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *vel_n, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *klo_vel_gb, const  int *khi_vel_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS</name>
      <anchor>a10</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb, const  int *klo_lse_rhs_gb, const  int *khi_lse_rhs_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *vel_n, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_ADD_CONST_CURV_TERM_TO_LSE_RHS</name>
      <anchor>a11</anchor>
      <arglist>(const  LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb, const  int *klo_lse_rhs_gb, const  int *khi_lse_rhs_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *phi_xx, const  LSMLIB_REAL *phi_xy, const  LSMLIB_REAL *phi_xz, const  LSMLIB_REAL *phi_yy, const  LSMLIB_REAL *phi_yz, const  LSMLIB_REAL *phi_zz, const  int *ilo_grad2_phi_gb, const  int *ihi_grad2_phi_gb, const  int *jlo_grad2_phi_gb, const  int *jhi_grad2_phi_gb, const  int *klo_grad2_phi_gb, const  int *khi_grad2_phi_gb, const  LSMLIB_REAL *b, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_ADD_CONST_PRECOMPUTED_CURV_TERM_TO_LSE_RHS</name>
      <anchor>a12</anchor>
      <arglist>(const  LSMLIB_REAL *lse_rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  int *klo_kappa_gb, const  int *khi_kappa_gb, const  LSMLIB_REAL *grad_mag_phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *b, const  int *ilo_rhs_fb, const  int *ihi_rhs_fb, const  int *jlo_rhs_fb, const  int *jhi_rhs_fb, const  int *klo_rhs_fb, const  int *khi_rhs_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_ADD_EXTERNAL_AND_NORMAL_VEL_TERM_TO_LSE_RHS</name>
      <anchor>a13</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *vel_n, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  LSMLIB_REAL *vel_z, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *klo_vel_gb, const  int *khi_vel_gb, const  int *ilo_rhs_fb, const  int *ihi_rhs_fb, const  int *jlo_rhs_fb, const  int *jhi_rhs_fb, const  int *klo_rhs_fb, const  int *khi_rhs_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_level_set_evolution3d_local.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/level_set_evolution/</path>
    <filename>lsm__level__set__evolution3d__local_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_ZERO_OUT_LEVEL_SET_EQN_RHS_LOCAL</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_ADD_ADVECTION_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_ADD_CONST_CURV_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_ADD_CONST_PRECOMPUTED_CURV_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_ADD_EXTERNAL_AND_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_ZERO_OUT_LEVEL_SET_EQN_RHS_LOCAL</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb, const  int *klo_lse_rhs_gb, const  int *khi_lse_rhs_gb, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_ADD_ADVECTION_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a8</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb, const  int *klo_lse_rhs_gb, const  int *khi_lse_rhs_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  LSMLIB_REAL *vel_z, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *klo_vel_gb, const  int *khi_vel_gb, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a9</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb, const  int *klo_lse_rhs_gb, const  int *khi_lse_rhs_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *vel_n, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *klo_vel_gb, const  int *khi_vel_gb, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_ADD_CONST_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a10</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb, const  int *klo_lse_rhs_gb, const  int *khi_lse_rhs_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *vel_n, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_ADD_CONST_CURV_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a11</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_lse_rhs_gb, const  int *ihi_lse_rhs_gb, const  int *jlo_lse_rhs_gb, const  int *jhi_lse_rhs_gb, const  int *klo_lse_rhs_gb, const  int *khi_lse_rhs_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *phi_xx, const  LSMLIB_REAL *phi_xy, const  LSMLIB_REAL *phi_xz, const  LSMLIB_REAL *phi_yy, const  LSMLIB_REAL *phi_yz, const  LSMLIB_REAL *phi_zz, const  int *ilo_grad2_phi_gb, const  int *ihi_grad2_phi_gb, const  int *jlo_grad2_phi_gb, const  int *jhi_grad2_phi_gb, const  int *klo_grad2_phi_gb, const  int *khi_grad2_phi_gb, const  LSMLIB_REAL *b, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_ADD_CONST_PRECOMPUTED_CURV_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a12</anchor>
      <arglist>(const  LSMLIB_REAL *lse_rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, LSMLIB_REAL *kappa, const  int *ilo_kappa_gb, const  int *ihi_kappa_gb, const  int *jlo_kappa_gb, const  int *jhi_kappa_gb, const  int *klo_kappa_gb, const  int *khi_kappa_gb, const  LSMLIB_REAL *grad_mag_phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *b, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_ADD_EXTERNAL_AND_NORMAL_VEL_TERM_TO_LSE_RHS_LOCAL</name>
      <anchor>a13</anchor>
      <arglist>(LSMLIB_REAL *lse_rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *vel_n, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  LSMLIB_REAL *vel_z, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *klo_vel_gb, const  int *khi_vel_gb, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_macros.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/serial/</path>
    <filename>lsm__macros_8h</filename>
    <includes id="lsm__grid_8h" name="lsm_grid.h" local="yes">lsm_grid.h</includes>
    <member kind="define">
      <type>#define</type>
      <name>SET_DATA_TO_CONSTANT</name>
      <anchor>a0</anchor>
      <arglist>(data, grid, value)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>NEGATE_DATA</name>
      <anchor>a1</anchor>
      <arglist>(data, grid)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IMPOSE_MASK</name>
      <anchor>a2</anchor>
      <arglist>(phi_masked, mask, phi, grid)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IMPOSE_MASK_LOCAL</name>
      <anchor>a3</anchor>
      <arglist>(phi_masked, mask, phi, grid, p)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IMPOSE_MIN</name>
      <anchor>a4</anchor>
      <arglist>(phi_min, phi1, phi2, grid)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>COPY_DATA</name>
      <anchor>a5</anchor>
      <arglist>(data_dst, data_src, grid)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>COMPUTE_MAX_ABS_ERR</name>
      <anchor>a6</anchor>
      <arglist>(max_abs_err, data1, data2, grid)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>COMPUTE_MAX_ABS_DATA</name>
      <anchor>a7</anchor>
      <arglist>(max_abs, data, grid)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>EXAMINE_ARRAY</name>
      <anchor>a8</anchor>
      <arglist>(name, data, g)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>PRINT_ARRAY_2D</name>
      <anchor>a9</anchor>
      <arglist>(data, grid)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_reinitialization1d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/reinitialization/</path>
    <filename>lsm__reinitialization1d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_COMPUTE_REINITIALIZATION_EQN_RHS</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_COMPUTE_REINITIALIZATION_EQN_RHS</name>
      <anchor>a1</anchor>
      <arglist>(LSMLIB_REAL *reinit_rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  LSMLIB_REAL *phi0, const  int *ilo_phi0_gb, const  int *ihi_phi0_gb, const  LSMLIB_REAL *phi_x_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dx, const  int *use_phi0_for_sgn)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_reinitialization2d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/reinitialization/</path>
    <filename>lsm__reinitialization2d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_REINITIALIZATION_EQN_RHS</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_ORTHOGONALIZATION_EQN_RHS</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_DISTANCE_FOR_SUBCELL_FIX</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_REINITIALIZATION_EQN_RHS_SUBCELL_FIX_ORDER1</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_REINITIALIZATION_EQN_RHS</name>
      <anchor>a4</anchor>
      <arglist>(LSMLIB_REAL *reinit_rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *phi0, const  int *ilo_phi0_gb, const  int *ihi_phi0_gb, const  int *jlo_phi0_gb, const  int *jhi_phi0_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *use_phi0_for_sgn)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_ORTHOGONALIZATION_EQN_RHS</name>
      <anchor>a5</anchor>
      <arglist>(LSMLIB_REAL *ortho_rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *psi, const  int *ilo_psi_gb, const  int *ihi_psi_gb, const  int *jlo_psi_gb, const  int *jhi_psi_gb, const  LSMLIB_REAL *psi_x_plus, const  LSMLIB_REAL *psi_y_plus, const  int *ilo_grad_psi_plus_gb, const  int *ihi_grad_psi_plus_gb, const  int *jlo_grad_psi_plus_gb, const  int *jhi_grad_psi_plus_gb, const  LSMLIB_REAL *psi_x_minus, const  LSMLIB_REAL *psi_y_minus, const  int *ilo_grad_psi_minus_gb, const  int *ihi_grad_psi_minus_gb, const  int *jlo_grad_psi_minus_gb, const  int *jhi_grad_psi_minus_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_DISTANCE_FOR_SUBCELL_FIX</name>
      <anchor>a6</anchor>
      <arglist>(LSMLIB_REAL *distance0, const  LSMLIB_REAL *phi0, const  int *ilo_phi0_gb, const  int *ihi_phi0_gb, const  int *jlo_phi0_gb, const  int *jhi_phi0_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_REINITIALIZATION_EQN_RHS_SUBCELL_FIX_ORDER1</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *reinit_rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *phi0, const  LSMLIB_REAL *distance0, const  int *ilo_phi0_gb, const  int *ihi_phi0_gb, const  int *jlo_phi0_gb, const  int *jhi_phi0_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_reinitialization2d_local.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/reinitialization/</path>
    <filename>lsm__reinitialization2d__local_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_REINITIALIZATION_EQN_RHS_LOCAL</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_ORTHOGONALIZATION_EQN_RHS_LOCAL</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_REINITIALIZATION_EQN_RHS_LOCAL</name>
      <anchor>a2</anchor>
      <arglist>(LSMLIB_REAL *reinit_rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *phi0, const  int *ilo_phi0_gb, const  int *ihi_phi0_gb, const  int *jlo_phi0_gb, const  int *jhi_phi0_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *use_phi0_for_sgn, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_ORTHOGONALIZATION_EQN_RHS_LOCAL</name>
      <anchor>a3</anchor>
      <arglist>(LSMLIB_REAL *ortho_rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *psi, const  int *ilo_psi_gb, const  int *ihi_psi_gb, const  int *jlo_psi_gb, const  int *jhi_psi_gb, const  LSMLIB_REAL *psi_x_plus, const  LSMLIB_REAL *psi_y_plus, const  int *ilo_grad_psi_plus_gb, const  int *ihi_grad_psi_plus_gb, const  int *jlo_grad_psi_plus_gb, const  int *jhi_grad_psi_plus_gb, const  LSMLIB_REAL *psi_x_minus, const  LSMLIB_REAL *psi_y_minus, const  int *ilo_grad_psi_minus_gb, const  int *ihi_grad_psi_minus_gb, const  int *jlo_grad_psi_minus_gb, const  int *jhi_grad_psi_minus_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_reinitialization3d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/reinitialization/</path>
    <filename>lsm__reinitialization3d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_REINITIALIZATION_EQN_RHS</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_ORTHOGONALIZATION_EQN_RHS</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_DISTANCE_FOR_SUBCELL_FIX</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_REINITIALIZATION_EQN_RHS_SUBCELL_FIX_ORDER1</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_REINITIALIZATION_EQN_RHS</name>
      <anchor>a4</anchor>
      <arglist>(LSMLIB_REAL *reinit_rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *phi0, const  int *ilo_phi0_gb, const  int *ihi_phi0_gb, const  int *jlo_phi0_gb, const  int *jhi_phi0_gb, const  int *klo_phi0_gb, const  int *khi_phi0_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  int *use_phi0_for_sgn)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_ORTHOGONALIZATION_EQN_RHS</name>
      <anchor>a5</anchor>
      <arglist>(LSMLIB_REAL *ortho_rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *psi, const  int *ilo_psi_gb, const  int *ihi_psi_gb, const  int *jlo_psi_gb, const  int *jhi_psi_gb, const  int *klo_psi_gb, const  int *khi_psi_gb, const  LSMLIB_REAL *psi_x_plus, const  LSMLIB_REAL *psi_y_plus, const  LSMLIB_REAL *psi_z_plus, const  int *ilo_grad_psi_plus_gb, const  int *ihi_grad_psi_plus_gb, const  int *jlo_grad_psi_plus_gb, const  int *jhi_grad_psi_plus_gb, const  int *klo_grad_psi_plus_gb, const  int *khi_grad_psi_plus_gb, const  LSMLIB_REAL *psi_x_minus, const  LSMLIB_REAL *psi_y_minus, const  LSMLIB_REAL *psi_z_minus, const  int *ilo_grad_psi_minus_gb, const  int *ihi_grad_psi_minus_gb, const  int *jlo_grad_psi_minus_gb, const  int *jhi_grad_psi_minus_gb, const  int *klo_grad_psi_minus_gb, const  int *khi_grad_psi_minus_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_DISTANCE_FOR_SUBCELL_FIX</name>
      <anchor>a6</anchor>
      <arglist>(LSMLIB_REAL *distance0, const  LSMLIB_REAL *phi0, const  int *ilo_phi0_gb, const  int *ihi_phi0_gb, const  int *jlo_phi0_gb, const  int *jhi_phi0_gb, const  int *klo_phi0_gb, const  int *khi_phi0_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_REINITIALIZATION_EQN_RHS_SUBCELL_FIX_ORDER1</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *reinit_rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *phi0, const  LSMLIB_REAL *distance0, const  int *ilo_phi0_gb, const  int *ihi_phi0_gb, const  int *jlo_phi0_gb, const  int *jhi_phi0_gb, const  int *klo_phi0_gb, const  int *khi_phi0_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_reinitialization3d_local.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/reinitialization/</path>
    <filename>lsm__reinitialization3d__local_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_REINITIALIZATION_EQN_RHS_LOCAL</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_ORTHOGONALIZATION_EQN_RHS_LOCAL</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_REINITIALIZATION_EQN_RHS_LOCAL</name>
      <anchor>a2</anchor>
      <arglist>(LSMLIB_REAL *reinit_rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *phi0, const  int *ilo_phi0_gb, const  int *ihi_phi0_gb, const  int *jlo_phi0_gb, const  int *jhi_phi0_gb, const  int *klo_phi0_gb, const  int *khi_phi0_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  int *use_phi0_for_sgn, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_ORTHOGONALIZATION_EQN_RHS_LOCAL</name>
      <anchor>a3</anchor>
      <arglist>(LSMLIB_REAL *ortho_rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *psi, const  int *ilo_psi_gb, const  int *ihi_psi_gb, const  int *jlo_psi_gb, const  int *jhi_psi_gb, const  int *klo_psi_gb, const  int *khi_psi_gb, const  LSMLIB_REAL *psi_x_plus, const  LSMLIB_REAL *psi_y_plus, const  LSMLIB_REAL *psi_z_plus, const  int *ilo_grad_psi_plus_gb, const  int *ihi_grad_psi_plus_gb, const  int *jlo_grad_psi_plus_gb, const  int *jhi_grad_psi_plus_gb, const  int *klo_grad_psi_plus_gb, const  int *khi_grad_psi_plus_gb, const  LSMLIB_REAL *psi_x_minus, const  LSMLIB_REAL *psi_y_minus, const  LSMLIB_REAL *psi_z_minus, const  int *ilo_grad_psi_minus_gb, const  int *ihi_grad_psi_minus_gb, const  int *jlo_grad_psi_minus_gb, const  int *jhi_grad_psi_minus_gb, const  int *klo_grad_psi_minus_gb, const  int *khi_grad_psi_minus_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_spatial_derivatives1d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/spatial_derivatives/</path>
    <filename>lsm__spatial__derivatives1d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_HJ_ENO1</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_HJ_ENO2</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_HJ_ENO3</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_HJ_WENO5</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_UPWIND_HJ_ENO1</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_UPWIND_HJ_ENO2</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_UPWIND_HJ_ENO3</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_UPWIND_HJ_WENO5</name>
      <anchor>a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_CENTRAL_GRAD_ORDER2</name>
      <anchor>a8</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_CENTRAL_GRAD_ORDER4</name>
      <anchor>a9</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_LAPLACIAN_ORDER2</name>
      <anchor>a10</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_PHI_UPWIND_GRAD_F</name>
      <anchor>a11</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_AVERAGE_GRAD_PHI</name>
      <anchor>a12</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_HJ_ENO1</name>
      <anchor>a13</anchor>
      <arglist>(LSMLIB_REAL *phi_x_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, LSMLIB_REAL *D1_x, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_HJ_ENO2</name>
      <anchor>a14</anchor>
      <arglist>(LSMLIB_REAL *phi_x_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, LSMLIB_REAL *D2, const  int *ilo_D2_gb, const  int *ihi_D2_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_HJ_ENO3</name>
      <anchor>a15</anchor>
      <arglist>(LSMLIB_REAL *phi_x_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, LSMLIB_REAL *D2, const  int *ilo_D2_gb, const  int *ihi_D2_gb, LSMLIB_REAL *D3, const  int *ilo_D3_gb, const  int *ihi_D3_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_HJ_WENO5</name>
      <anchor>a16</anchor>
      <arglist>(LSMLIB_REAL *phi_x_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_UPWIND_HJ_ENO1</name>
      <anchor>a17</anchor>
      <arglist>(LSMLIB_REAL *phi_x, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  LSMLIB_REAL *vel_x, const  int *ilo_vel_gb, const  int *ihi_vel_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_UPWIND_HJ_ENO2</name>
      <anchor>a18</anchor>
      <arglist>(LSMLIB_REAL *phi_x, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  LSMLIB_REAL *vel_x, const  int *ilo_vel_gb, const  int *ihi_vel_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, LSMLIB_REAL *D2, const  int *ilo_D2_gb, const  int *ihi_D2_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_UPWIND_HJ_ENO3</name>
      <anchor>a19</anchor>
      <arglist>(LSMLIB_REAL *phi_x, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  LSMLIB_REAL *vel_x, const  int *ilo_vel_gb, const  int *ihi_vel_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, LSMLIB_REAL *D2, const  int *ilo_D2_gb, const  int *ihi_D2_gb, LSMLIB_REAL *D3, const  int *ilo_D3_gb, const  int *ihi_D3_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_UPWIND_HJ_WENO5</name>
      <anchor>a20</anchor>
      <arglist>(LSMLIB_REAL *phi_x, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  LSMLIB_REAL *vel_x, const  int *ilo_vel_gb, const  int *ihi_vel_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_CENTRAL_GRAD_ORDER2</name>
      <anchor>a21</anchor>
      <arglist>(LSMLIB_REAL *phi_x, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_CENTRAL_GRAD_ORDER4</name>
      <anchor>a22</anchor>
      <arglist>(LSMLIB_REAL *phi_x, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_LAPLACIAN_ORDER2</name>
      <anchor>a23</anchor>
      <arglist>(LSMLIB_REAL *laplacian_phi, const  int *ilo_laplacian_phi_gb, const  int *ihi_laplacian_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dx)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_PHI_UPWIND_GRAD_F</name>
      <anchor>a24</anchor>
      <arglist>(LSMLIB_REAL *F_x, const  int *ilo_grad_F_gb, const  int *ihi_grad_F_gb, LSMLIB_REAL *F_x_plus, const  int *ilo_grad_F_plus_gb, const  int *ihi_grad_F_plus_gb, LSMLIB_REAL *F_x_minus, const  int *ilo_grad_F_minus_gb, const  int *ihi_grad_F_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *ilo_fb, const  int *ihi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_AVERAGE_GRAD_PHI</name>
      <anchor>a25</anchor>
      <arglist>(LSMLIB_REAL *phi_x, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, LSMLIB_REAL *phi_x_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *ilo_fb, const  int *ihi_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_spatial_derivatives2d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/spatial_derivatives/</path>
    <filename>lsm__spatial__derivatives2d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_HJ_ENO1</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_HJ_ENO2</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_HJ_ENO3</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_HJ_WENO5</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_UPWIND_HJ_ENO1</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_UPWIND_HJ_ENO2</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_UPWIND_HJ_ENO3</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_UPWIND_HJ_WENO5</name>
      <anchor>a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_CENTRAL_GRAD_ORDER2</name>
      <anchor>a8</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_CENTRAL_GRAD_ORDER4</name>
      <anchor>a9</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_LAPLACIAN_ORDER2</name>
      <anchor>a10</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_PHI_UPWIND_GRAD_F</name>
      <anchor>a11</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_AVERAGE_GRAD_PHI</name>
      <anchor>a12</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_GRADIENT_MAGNITUDE</name>
      <anchor>a13</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_DIVERGENCE_CENTRAL</name>
      <anchor>a14</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_HJ_ENO1</name>
      <anchor>a15</anchor>
      <arglist>(LSMLIB_REAL *phi_x_plus, LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_HJ_ENO2</name>
      <anchor>a16</anchor>
      <arglist>(LSMLIB_REAL *phi_x_plus, LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, LSMLIB_REAL *D2, const  int *ilo_D2_gb, const  int *ihi_D2_gb, const  int *jlo_D2_gb, const  int *jhi_D2_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_HJ_ENO3</name>
      <anchor>a17</anchor>
      <arglist>(LSMLIB_REAL *phi_x_plus, LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, LSMLIB_REAL *D2, const  int *ilo_D2_gb, const  int *ihi_D2_gb, const  int *jlo_D2_gb, const  int *jhi_D2_gb, LSMLIB_REAL *D3, const  int *ilo_D3_gb, const  int *ihi_D3_gb, const  int *jlo_D3_gb, const  int *jhi_D3_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_HJ_WENO5</name>
      <anchor>a18</anchor>
      <arglist>(LSMLIB_REAL *phi_x_plus, LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_UPWIND_HJ_ENO1</name>
      <anchor>a19</anchor>
      <arglist>(LSMLIB_REAL *phi_x, LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_UPWIND_HJ_ENO2</name>
      <anchor>a20</anchor>
      <arglist>(LSMLIB_REAL *phi_x, LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, LSMLIB_REAL *D2, const  int *ilo_D2_gb, const  int *ihi_D2_gb, const  int *jlo_D2_gb, const  int *jhi_D2_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_UPWIND_HJ_ENO3</name>
      <anchor>a21</anchor>
      <arglist>(LSMLIB_REAL *phi_x, LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, LSMLIB_REAL *D2, const  int *ilo_D2_gb, const  int *ihi_D2_gb, const  int *jlo_D2_gb, const  int *jhi_D2_gb, LSMLIB_REAL *D3, const  int *ilo_D3_gb, const  int *ihi_D3_gb, const  int *jlo_D3_gb, const  int *jhi_D3_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_UPWIND_HJ_WENO5</name>
      <anchor>a22</anchor>
      <arglist>(LSMLIB_REAL *phi_x, LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_CENTRAL_GRAD_ORDER2</name>
      <anchor>a23</anchor>
      <arglist>(LSMLIB_REAL *phi_x, LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_CENTRAL_GRAD_ORDER4</name>
      <anchor>a24</anchor>
      <arglist>(LSMLIB_REAL *phi_x, LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_LAPLACIAN_ORDER2</name>
      <anchor>a25</anchor>
      <arglist>(LSMLIB_REAL *laplacian_phi, const  int *ilo_laplacian_phi_gb, const  int *ihi_laplacian_phi_gb, const  int *jlo_laplacian_phi_gb, const  int *jhi_laplacian_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_PHI_UPWIND_GRAD_F</name>
      <anchor>a26</anchor>
      <arglist>(LSMLIB_REAL *F_x, LSMLIB_REAL *F_y, const  int *ilo_grad_F_gb, const  int *ihi_grad_F_gb, const  int *jlo_grad_F_gb, const  int *jhi_grad_F_gb, LSMLIB_REAL *F_x_plus, LSMLIB_REAL *F_y_plus, const  int *ilo_grad_F_plus_gb, const  int *ihi_grad_F_plus_gb, const  int *jlo_grad_F_plus_gb, const  int *jhi_grad_F_plus_gb, LSMLIB_REAL *F_x_minus, LSMLIB_REAL *F_y_minus, const  int *ilo_grad_F_minus_gb, const  int *ihi_grad_F_minus_gb, const  int *jlo_grad_F_minus_gb, const  int *jhi_grad_F_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_AVERAGE_GRAD_PHI</name>
      <anchor>a27</anchor>
      <arglist>(LSMLIB_REAL *phi_x, LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, LSMLIB_REAL *phi_x_plus, LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_GRADIENT_MAGNITUDE</name>
      <anchor>a28</anchor>
      <arglist>(const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_DIVERGENCE_CENTRAL</name>
      <anchor>a29</anchor>
      <arglist>(LSMLIB_REAL *divF, const  int *ilo_divf_gb, const  int *ihi_divf_gb, const  int *jlo_divf_gb, const  int *jhi_divf_gb, const  LSMLIB_REAL *FX, const  LSMLIB_REAL *FY, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_spatial_derivatives2d_local.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/spatial_derivatives/</path>
    <filename>lsm__spatial__derivatives2d__local_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_HJ_ENO1_LOCAL</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_HJ_ENO2_LOCAL</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_HJ_ENO3_LOCAL</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_HJ_WENO5_LOCAL</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_UPWIND_HJ_ENO2_LOCAL</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_CENTRAL_GRAD_ORDER2_LOCAL</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_CENTRAL_GRAD_ORDER4_LOCAL</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_LAPLACIAN_ORDER2_LOCAL</name>
      <anchor>a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_AVE_GRAD_PHI_LOCAL</name>
      <anchor>a8</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_GRADIENT_MAGNITUDE_LOCAL</name>
      <anchor>a9</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_DIVERGENCE_CENTRAL_LOCAL</name>
      <anchor>a10</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_HJ_ENO1_LOCAL</name>
      <anchor>a11</anchor>
      <arglist>(LSMLIB_REAL *phi_x_plus, LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index0, const  int *nhi_index0, const  int *nlo_index1, const  int *nhi_index1, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb, const  unsigned char *mark_D1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_HJ_ENO2_LOCAL</name>
      <anchor>a12</anchor>
      <arglist>(LSMLIB_REAL *phi_x_plus, LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, LSMLIB_REAL *D2, const  int *ilo_D2_gb, const  int *ihi_D2_gb, const  int *jlo_D2_gb, const  int *jhi_D2_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index0, const  int *nhi_index0, const  int *nlo_index1, const  int *nhi_index1, const  int *nlo_index2, const  int *nhi_index2, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb, const  unsigned char *mark_D1, const  unsigned char *mark_D2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_HJ_ENO3_LOCAL</name>
      <anchor>a13</anchor>
      <arglist>(LSMLIB_REAL *phi_x_plus, LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, LSMLIB_REAL *D2, const  int *ilo_D2_gb, const  int *ihi_D2_gb, const  int *jlo_D2_gb, const  int *jhi_D2_gb, LSMLIB_REAL *D3, const  int *ilo_D3_gb, const  int *ihi_D3_gb, const  int *jlo_D3_gb, const  int *jhi_D3_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index0, const  int *nhi_index0, const  int *nlo_index1, const  int *nhi_index1, const  int *nlo_index2, const  int *nhi_index2, const  int *nlo_index3, const  int *nhi_index3, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb, const  unsigned char *mark_D1, const  unsigned char *mark_D2, const  unsigned char *mark_D3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_HJ_WENO5_LOCAL</name>
      <anchor>a14</anchor>
      <arglist>(LSMLIB_REAL *phi_x_plus, LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index0, const  int *nhi_index0, const  int *nlo_index1, const  int *nhi_index1, const  int *nlo_index2, const  int *nhi_index2, const  int *nlo_index3, const  int *nhi_index3, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb, const  unsigned char *mark_D1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_UPWIND_HJ_ENO2_LOCAL</name>
      <anchor>a15</anchor>
      <arglist>(LSMLIB_REAL *phi_x, LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, LSMLIB_REAL *D2, const  int *ilo_D2_gb, const  int *ihi_D2_gb, const  int *jlo_D2_gb, const  int *jhi_D2_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index0, const  int *nhi_index0, const  int *nlo_index1, const  int *nhi_index1, const  int *nlo_index2, const  int *nhi_index2, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb, const  unsigned char *mark_D1, const  unsigned char *mark_D2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_CENTRAL_GRAD_ORDER2_LOCAL</name>
      <anchor>a16</anchor>
      <arglist>(LSMLIB_REAL *phi_x, LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_CENTRAL_GRAD_ORDER4_LOCAL</name>
      <anchor>a17</anchor>
      <arglist>(LSMLIB_REAL *phi_x, LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_LAPLACIAN_ORDER2_LOCAL</name>
      <anchor>a18</anchor>
      <arglist>(LSMLIB_REAL *laplacian_phi, const  int *ilo_laplacian_phi_gb, const  int *ihi_laplacian_phi_gb, const  int *jlo_laplacian_phi_gb, const  int *jhi_laplacian_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_AVE_GRAD_PHI_LOCAL</name>
      <anchor>a19</anchor>
      <arglist>(LSMLIB_REAL *grad_phi_ave, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_GRADIENT_MAGNITUDE_LOCAL</name>
      <anchor>a20</anchor>
      <arglist>(const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_DIVERGENCE_CENTRAL_LOCAL</name>
      <anchor>a21</anchor>
      <arglist>(LSMLIB_REAL *divF, const  int *ilo_divf_gb, const  int *ihi_divf_gb, const  int *jlo_divf_gb, const  int *jhi_divf_gb, const  LSMLIB_REAL *FX, const  LSMLIB_REAL *FY, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_spatial_derivatives3d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/spatial_derivatives/</path>
    <filename>lsm__spatial__derivatives3d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_HJ_ENO1</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_HJ_ENO2</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_HJ_ENO3</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_HJ_WENO5</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_UPWIND_HJ_ENO1</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_UPWIND_HJ_ENO2</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_UPWIND_HJ_ENO3</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_UPWIND_HJ_WENO5</name>
      <anchor>a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_CENTRAL_GRAD_ORDER2</name>
      <anchor>a8</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_CENTRAL_GRAD_ORDER4</name>
      <anchor>a9</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_LAPLACIAN_ORDER2</name>
      <anchor>a10</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_PHI_UPWIND_GRAD_F</name>
      <anchor>a11</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_AVERAGE_GRAD_PHI</name>
      <anchor>a12</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_GRADIENT_MAGNITUDE</name>
      <anchor>a13</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_HJ_ENO1</name>
      <anchor>a14</anchor>
      <arglist>(LSMLIB_REAL *phi_x_plus, LSMLIB_REAL *phi_y_plus, LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, LSMLIB_REAL *phi_y_minus, LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, const  int *klo_D1_gb, const  int *khi_D1_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_HJ_ENO2</name>
      <anchor>a15</anchor>
      <arglist>(LSMLIB_REAL *phi_x_plus, LSMLIB_REAL *phi_y_plus, LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, LSMLIB_REAL *phi_y_minus, LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, const  int *klo_D1_gb, const  int *khi_D1_gb, LSMLIB_REAL *D2, const  int *ilo_D2_gb, const  int *ihi_D2_gb, const  int *jlo_D2_gb, const  int *jhi_D2_gb, const  int *klo_D2_gb, const  int *khi_D2_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_HJ_ENO3</name>
      <anchor>a16</anchor>
      <arglist>(LSMLIB_REAL *phi_x_plus, LSMLIB_REAL *phi_y_plus, LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, LSMLIB_REAL *phi_y_minus, LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, const  int *klo_D1_gb, const  int *khi_D1_gb, LSMLIB_REAL *D2, const  int *ilo_D2_gb, const  int *ihi_D2_gb, const  int *jlo_D2_gb, const  int *jhi_D2_gb, const  int *klo_D2_gb, const  int *khi_D2_gb, LSMLIB_REAL *D3, const  int *ilo_D3_gb, const  int *ihi_D3_gb, const  int *jlo_D3_gb, const  int *jhi_D3_gb, const  int *klo_D3_gb, const  int *khi_D3_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_HJ_WENO5</name>
      <anchor>a17</anchor>
      <arglist>(LSMLIB_REAL *phi_x_plus, LSMLIB_REAL *phi_y_plus, LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, LSMLIB_REAL *phi_y_minus, LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, const  int *klo_D1_gb, const  int *khi_D1_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_UPWIND_HJ_ENO1</name>
      <anchor>a18</anchor>
      <arglist>(LSMLIB_REAL *phi_x, LSMLIB_REAL *phi_y, LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  LSMLIB_REAL *vel_z, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *klo_vel_gb, const  int *khi_vel_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, const  int *klo_D1_gb, const  int *khi_D1_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_UPWIND_HJ_ENO2</name>
      <anchor>a19</anchor>
      <arglist>(LSMLIB_REAL *phi_x, LSMLIB_REAL *phi_y, LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  LSMLIB_REAL *vel_z, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *klo_vel_gb, const  int *khi_vel_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, const  int *klo_D1_gb, const  int *khi_D1_gb, LSMLIB_REAL *D2, const  int *ilo_D2_gb, const  int *ihi_D2_gb, const  int *jlo_D2_gb, const  int *jhi_D2_gb, const  int *klo_D2_gb, const  int *khi_D2_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_UPWIND_HJ_ENO3</name>
      <anchor>a20</anchor>
      <arglist>(LSMLIB_REAL *phi_x, LSMLIB_REAL *phi_y, LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  LSMLIB_REAL *vel_z, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *klo_vel_gb, const  int *khi_vel_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, const  int *klo_D1_gb, const  int *khi_D1_gb, LSMLIB_REAL *D2, const  int *ilo_D2_gb, const  int *ihi_D2_gb, const  int *jlo_D2_gb, const  int *jhi_D2_gb, const  int *klo_D2_gb, const  int *khi_D2_gb, LSMLIB_REAL *D3, const  int *ilo_D3_gb, const  int *ihi_D3_gb, const  int *jlo_D3_gb, const  int *jhi_D3_gb, const  int *klo_D3_gb, const  int *khi_D3_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_UPWIND_HJ_WENO5</name>
      <anchor>a21</anchor>
      <arglist>(LSMLIB_REAL *phi_x, LSMLIB_REAL *phi_y, LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  LSMLIB_REAL *vel_z, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *klo_vel_gb, const  int *khi_vel_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, const  int *klo_D1_gb, const  int *khi_D1_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_CENTRAL_GRAD_ORDER2</name>
      <anchor>a22</anchor>
      <arglist>(LSMLIB_REAL *phi_x, LSMLIB_REAL *phi_y, LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_CENTRAL_GRAD_ORDER4</name>
      <anchor>a23</anchor>
      <arglist>(LSMLIB_REAL *phi_x, LSMLIB_REAL *phi_y, LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_LAPLACIAN_ORDER2</name>
      <anchor>a24</anchor>
      <arglist>(LSMLIB_REAL *laplacian_phi, const  int *ilo_laplacian_phi_gb, const  int *ihi_laplacian_phi_gb, const  int *jlo_laplacian_phi_gb, const  int *jhi_laplacian_phi_gb, const  int *klo_laplacian_phi_gb, const  int *khi_laplacian_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_PHI_UPWIND_GRAD_F</name>
      <anchor>a25</anchor>
      <arglist>(LSMLIB_REAL *F_x, LSMLIB_REAL *F_y, LSMLIB_REAL *F_z, const  int *ilo_grad_F_gb, const  int *ihi_grad_F_gb, const  int *jlo_grad_F_gb, const  int *jhi_grad_F_gb, const  int *klo_grad_F_gb, const  int *khi_grad_F_gb, LSMLIB_REAL *F_x_plus, LSMLIB_REAL *F_y_plus, LSMLIB_REAL *F_z_plus, const  int *ilo_grad_F_plus_gb, const  int *ihi_grad_F_plus_gb, const  int *jlo_grad_F_plus_gb, const  int *jhi_grad_F_plus_gb, const  int *klo_grad_F_plus_gb, const  int *khi_grad_F_plus_gb, LSMLIB_REAL *F_x_minus, LSMLIB_REAL *F_y_minus, LSMLIB_REAL *F_z_minus, const  int *ilo_grad_F_minus_gb, const  int *ihi_grad_F_minus_gb, const  int *jlo_grad_F_minus_gb, const  int *jhi_grad_F_minus_gb, const  int *klo_grad_F_minus_gb, const  int *khi_grad_F_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_AVERAGE_GRAD_PHI</name>
      <anchor>a26</anchor>
      <arglist>(LSMLIB_REAL *phi_x, LSMLIB_REAL *phi_y, LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, LSMLIB_REAL *phi_x_plus, LSMLIB_REAL *phi_y_plus, LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, LSMLIB_REAL *phi_y_minus, LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_GRADIENT_MAGNITUDE</name>
      <anchor>a27</anchor>
      <arglist>(const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  int *ilo_grad_phi_fb, const  int *ihi_grad_phi_fb, const  int *jlo_grad_phi_fb, const  int *jhi_grad_phi_fb, const  int *klo_grad_phi_fb, const  int *khi_grad_phi_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_spatial_derivatives3d_local.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/spatial_derivatives/</path>
    <filename>lsm__spatial__derivatives3d__local_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_HJ_ENO1_LOCAL</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_HJ_ENO2_LOCAL</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_CENTRAL_GRAD_ORDER2_LOCAL</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_CENTRAL_GRAD_ORDER4_LOCAL</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_LAPLACIAN_ORDER2_LOCAL</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_AVE_GRAD_PHI_LOCAL</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_GRADIENT_MAGNITUDE_LOCAL</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_HJ_ENO1_LOCAL</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *phi_x_plus, LSMLIB_REAL *phi_y_plus, LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, LSMLIB_REAL *phi_y_minus, LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, const  int *klo_D1_gb, const  int *khi_D1_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index0, const  int *nhi_index0, const  int *nlo_index1, const  int *nhi_index1, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb, const  unsigned char *mark_D1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_HJ_ENO2_LOCAL</name>
      <anchor>a8</anchor>
      <arglist>(LSMLIB_REAL *phi_x_plus, LSMLIB_REAL *phi_y_plus, LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, LSMLIB_REAL *phi_x_minus, LSMLIB_REAL *phi_y_minus, LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, LSMLIB_REAL *D1, const  int *ilo_D1_gb, const  int *ihi_D1_gb, const  int *jlo_D1_gb, const  int *jhi_D1_gb, const  int *klo_D1_gb, const  int *khi_D1_gb, LSMLIB_REAL *D2, const  int *ilo_D2_gb, const  int *ihi_D2_gb, const  int *jlo_D2_gb, const  int *jhi_D2_gb, const  int *klo_D2_gb, const  int *khi_D2_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index0, const  int *nhi_index0, const  int *nlo_index1, const  int *nhi_index1, const  int *nlo_index2, const  int *nhi_index2, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb, const  unsigned char *mark_D1, const  unsigned char *mark_D2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_CENTRAL_GRAD_ORDER2_LOCAL</name>
      <anchor>a9</anchor>
      <arglist>(LSMLIB_REAL *phi_x, LSMLIB_REAL *phi_y, LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_CENTRAL_GRAD_ORDER4_LOCAL</name>
      <anchor>a10</anchor>
      <arglist>(LSMLIB_REAL *phi_x, LSMLIB_REAL *phi_y, LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_LAPLACIAN_ORDER2_LOCAL</name>
      <anchor>a11</anchor>
      <arglist>(LSMLIB_REAL *laplacian_phi, const  int *ilo_laplacian_phi_gb, const  int *ihi_laplacian_phi_gb, const  int *jlo_laplacian_phi_gb, const  int *jhi_laplacian_phi_gb, const  int *klo_laplacian_phi_gb, const  int *khi_laplacian_phi_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_AVE_GRAD_PHI_LOCAL</name>
      <anchor>a12</anchor>
      <arglist>(LSMLIB_REAL *grad_phi_ave, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_GRADIENT_MAGNITUDE_LOCAL</name>
      <anchor>a13</anchor>
      <arglist>(const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_tvd_runge_kutta1d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/time_integration/</path>
    <filename>lsm__tvd__runge__kutta1d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_RK1_STEP</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_TVD_RK2_STAGE1</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_TVD_RK2_STAGE2</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_TVD_RK3_STAGE1</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_TVD_RK3_STAGE2</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_TVD_RK3_STAGE3</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_RK1_STEP</name>
      <anchor>a6</anchor>
      <arglist>(LSMLIB_REAL *u_next, const  int *ilo_u_next_gb, const  int *ihi_u_next_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_TVD_RK2_STAGE1</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_TVD_RK2_STAGE2</name>
      <anchor>a8</anchor>
      <arglist>(LSMLIB_REAL *u_next, const  int *ilo_u_next_gb, const  int *ihi_u_next_gb, const  LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_TVD_RK3_STAGE1</name>
      <anchor>a9</anchor>
      <arglist>(LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_TVD_RK3_STAGE2</name>
      <anchor>a10</anchor>
      <arglist>(LSMLIB_REAL *u_stage2, const  int *ilo_u_stage2_gb, const  int *ihi_u_stage2_gb, const  LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_TVD_RK3_STAGE3</name>
      <anchor>a11</anchor>
      <arglist>(LSMLIB_REAL *u_next, const  int *ilo_u_next_gb, const  int *ihi_u_next_gb, const  LSMLIB_REAL *u_stage2, const  int *ilo_u_stage2_gb, const  int *ihi_u_stage2_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *ilo_fb, const  int *ihi_fb, const  LSMLIB_REAL *dt)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_tvd_runge_kutta2d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/time_integration/</path>
    <filename>lsm__tvd__runge__kutta2d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_RK1_STEP</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_TVD_RK2_STAGE1</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_TVD_RK2_STAGE2</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_TVD_RK3_STAGE1</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_TVD_RK3_STAGE2</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_TVD_RK3_STAGE3</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_RK1_STEP</name>
      <anchor>a6</anchor>
      <arglist>(LSMLIB_REAL *u_next, const  int *ilo_u_next_gb, const  int *ihi_u_next_gb, const  int *jlo_u_next_gb, const  int *jhi_u_next_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_TVD_RK2_STAGE1</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  int *jlo_u_stage1_gb, const  int *jhi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_TVD_RK2_STAGE2</name>
      <anchor>a8</anchor>
      <arglist>(LSMLIB_REAL *u_next, const  int *ilo_u_next_gb, const  int *ihi_u_next_gb, const  int *jlo_u_next_gb, const  int *jhi_u_next_gb, const  LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  int *jlo_u_stage1_gb, const  int *jhi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_TVD_RK3_STAGE1</name>
      <anchor>a9</anchor>
      <arglist>(LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  int *jlo_u_stage1_gb, const  int *jhi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_TVD_RK3_STAGE2</name>
      <anchor>a10</anchor>
      <arglist>(LSMLIB_REAL *u_stage2, const  int *ilo_u_stage2_gb, const  int *ihi_u_stage2_gb, const  int *jlo_u_stage2_gb, const  int *jhi_u_stage2_gb, const  LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  int *jlo_u_stage1_gb, const  int *jhi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_TVD_RK3_STAGE3</name>
      <anchor>a11</anchor>
      <arglist>(LSMLIB_REAL *u_next, const  int *ilo_u_next_gb, const  int *ihi_u_next_gb, const  int *jlo_u_next_gb, const  int *jhi_u_next_gb, const  LSMLIB_REAL *u_stage2, const  int *ilo_u_stage2_gb, const  int *ihi_u_stage2_gb, const  int *jlo_u_stage2_gb, const  int *jhi_u_stage2_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  LSMLIB_REAL *dt)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_tvd_runge_kutta2d_local.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/time_integration/</path>
    <filename>lsm__tvd__runge__kutta2d__local_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_RK1_STEP_LOCAL</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_TVD_RK2_STAGE1_LOCAL</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_TVD_RK2_STAGE2_LOCAL</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_TVD_RK3_STAGE1_LOCAL</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_TVD_RK3_STAGE2_LOCAL</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_TVD_RK3_STAGE3_LOCAL</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_RK1_STEP_LOCAL</name>
      <anchor>a6</anchor>
      <arglist>(LSMLIB_REAL *u_next, const  int *ilo_u_next_gb, const  int *ihi_u_next_gb, const  int *jlo_u_next_gb, const  int *jhi_u_next_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  LSMLIB_REAL *dt, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_TVD_RK2_STAGE1_LOCAL</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  int *jlo_u_stage1_gb, const  int *jhi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  LSMLIB_REAL *dt, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_TVD_RK2_STAGE2_LOCAL</name>
      <anchor>a8</anchor>
      <arglist>(LSMLIB_REAL *u_next, const  int *ilo_u_next_gb, const  int *ihi_u_next_gb, const  int *jlo_u_next_gb, const  int *jhi_u_next_gb, const  LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  int *jlo_u_stage1_gb, const  int *jhi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  LSMLIB_REAL *dt, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_TVD_RK3_STAGE1_LOCAL</name>
      <anchor>a9</anchor>
      <arglist>(LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  int *jlo_u_stage1_gb, const  int *jhi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  LSMLIB_REAL *dt, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_TVD_RK3_STAGE2_LOCAL</name>
      <anchor>a10</anchor>
      <arglist>(LSMLIB_REAL *u_stage2, const  int *ilo_u_stage2_gb, const  int *ihi_u_stage2_gb, const  int *jlo_u_stage2_gb, const  int *jhi_u_stage2_gb, const  LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  int *jlo_u_stage1_gb, const  int *jhi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  LSMLIB_REAL *dt, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_TVD_RK3_STAGE3_LOCAL</name>
      <anchor>a11</anchor>
      <arglist>(LSMLIB_REAL *u_next, const  int *ilo_u_next_gb, const  int *ihi_u_next_gb, const  int *jlo_u_next_gb, const  int *jhi_u_next_gb, const  LSMLIB_REAL *u_stage2, const  int *ilo_u_stage2_gb, const  int *ihi_u_stage2_gb, const  int *jlo_u_stage2_gb, const  int *jhi_u_stage2_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  LSMLIB_REAL *dt, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_tvd_runge_kutta3d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/time_integration/</path>
    <filename>lsm__tvd__runge__kutta3d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_RK1_STEP</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_TVD_RK2_STAGE1</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_TVD_RK2_STAGE2</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_TVD_RK3_STAGE1</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_TVD_RK3_STAGE2</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_TVD_RK3_STAGE3</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_RK1_STEP</name>
      <anchor>a6</anchor>
      <arglist>(LSMLIB_REAL *u_next, const  int *ilo_u_next_gb, const  int *ihi_u_next_gb, const  int *jlo_u_next_gb, const  int *jhi_u_next_gb, const  int *klo_u_next_gb, const  int *khi_u_next_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  int *klo_u_cur_gb, const  int *khi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_TVD_RK2_STAGE1</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  int *jlo_u_stage1_gb, const  int *jhi_u_stage1_gb, const  int *klo_u_stage1_gb, const  int *khi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  int *klo_u_cur_gb, const  int *khi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_TVD_RK2_STAGE2</name>
      <anchor>a8</anchor>
      <arglist>(LSMLIB_REAL *u_next, const  int *ilo_u_next_gb, const  int *ihi_u_next_gb, const  int *jlo_u_next_gb, const  int *jhi_u_next_gb, const  int *klo_u_next_gb, const  int *khi_u_next_gb, const  LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  int *jlo_u_stage1_gb, const  int *jhi_u_stage1_gb, const  int *klo_u_stage1_gb, const  int *khi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  int *klo_u_cur_gb, const  int *khi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_TVD_RK3_STAGE1</name>
      <anchor>a9</anchor>
      <arglist>(LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  int *jlo_u_stage1_gb, const  int *jhi_u_stage1_gb, const  int *klo_u_stage1_gb, const  int *khi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  int *klo_u_cur_gb, const  int *khi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_TVD_RK3_STAGE2</name>
      <anchor>a10</anchor>
      <arglist>(LSMLIB_REAL *u_stage2, const  int *ilo_u_stage2_gb, const  int *ihi_u_stage2_gb, const  int *jlo_u_stage2_gb, const  int *jhi_u_stage2_gb, const  int *klo_u_stage2_gb, const  int *khi_u_stage2_gb, const  LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  int *jlo_u_stage1_gb, const  int *jhi_u_stage1_gb, const  int *klo_u_stage1_gb, const  int *khi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  int *klo_u_cur_gb, const  int *khi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_TVD_RK3_STAGE3</name>
      <anchor>a11</anchor>
      <arglist>(LSMLIB_REAL *u_next, const  int *ilo_u_next_gb, const  int *ihi_u_next_gb, const  int *jlo_u_next_gb, const  int *jhi_u_next_gb, const  int *klo_u_next_gb, const  int *khi_u_next_gb, const  LSMLIB_REAL *u_stage2, const  int *ilo_u_stage2_gb, const  int *ihi_u_stage2_gb, const  int *jlo_u_stage2_gb, const  int *jhi_u_stage2_gb, const  int *klo_u_stage2_gb, const  int *khi_u_stage2_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  int *klo_u_cur_gb, const  int *khi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb, const  LSMLIB_REAL *dt)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_tvd_runge_kutta3d_local.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/time_integration/</path>
    <filename>lsm__tvd__runge__kutta3d__local_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_RK1_STEP_LOCAL</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_TVD_RK2_STAGE1_LOCAL</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_TVD_RK2_STAGE2_LOCAL</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_TVD_RK3_STAGE1_LOCAL</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_TVD_RK3_STAGE2_LOCAL</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_TVD_RK3_STAGE3_LOCAL</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_RK1_STEP_LOCAL</name>
      <anchor>a6</anchor>
      <arglist>(LSMLIB_REAL *u_next, const  int *ilo_u_next_gb, const  int *ihi_u_next_gb, const  int *jlo_u_next_gb, const  int *jhi_u_next_gb, const  int *klo_u_next_gb, const  int *khi_u_next_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  int *klo_u_cur_gb, const  int *khi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  LSMLIB_REAL *dt, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_TVD_RK2_STAGE1_LOCAL</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  int *jlo_u_stage1_gb, const  int *jhi_u_stage1_gb, const  int *klo_u_stage1_gb, const  int *khi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  int *klo_u_cur_gb, const  int *khi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  LSMLIB_REAL *dt, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_TVD_RK2_STAGE2_LOCAL</name>
      <anchor>a8</anchor>
      <arglist>(LSMLIB_REAL *u_next, const  int *ilo_u_next_gb, const  int *ihi_u_next_gb, const  int *jlo_u_next_gb, const  int *jhi_u_next_gb, const  int *klo_u_next_gb, const  int *khi_u_next_gb, const  LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  int *jlo_u_stage1_gb, const  int *jhi_u_stage1_gb, const  int *klo_u_stage1_gb, const  int *khi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  int *klo_u_cur_gb, const  int *khi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  LSMLIB_REAL *dt, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_TVD_RK3_STAGE1_LOCAL</name>
      <anchor>a9</anchor>
      <arglist>(LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  int *jlo_u_stage1_gb, const  int *jhi_u_stage1_gb, const  int *klo_u_stage1_gb, const  int *khi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  int *klo_u_cur_gb, const  int *khi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  LSMLIB_REAL *dt, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_TVD_RK3_STAGE2_LOCAL</name>
      <anchor>a10</anchor>
      <arglist>(LSMLIB_REAL *u_stage2, const  int *ilo_u_stage2_gb, const  int *ihi_u_stage2_gb, const  int *jlo_u_stage2_gb, const  int *jhi_u_stage2_gb, const  int *klo_u_stage2_gb, const  int *khi_u_stage2_gb, const  LSMLIB_REAL *u_stage1, const  int *ilo_u_stage1_gb, const  int *ihi_u_stage1_gb, const  int *jlo_u_stage1_gb, const  int *jhi_u_stage1_gb, const  int *klo_u_stage1_gb, const  int *khi_u_stage1_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  int *klo_u_cur_gb, const  int *khi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  LSMLIB_REAL *dt, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_TVD_RK3_STAGE3_LOCAL</name>
      <anchor>a11</anchor>
      <arglist>(LSMLIB_REAL *u_next, const  int *ilo_u_next_gb, const  int *ihi_u_next_gb, const  int *jlo_u_next_gb, const  int *jhi_u_next_gb, const  int *klo_u_next_gb, const  int *khi_u_next_gb, const  LSMLIB_REAL *u_stage2, const  int *ilo_u_stage2_gb, const  int *ihi_u_stage2_gb, const  int *jlo_u_stage2_gb, const  int *jhi_u_stage2_gb, const  int *klo_u_stage2_gb, const  int *khi_u_stage2_gb, const  LSMLIB_REAL *u_cur, const  int *ilo_u_cur_gb, const  int *ihi_u_cur_gb, const  int *jlo_u_cur_gb, const  int *jhi_u_cur_gb, const  int *klo_u_cur_gb, const  int *khi_u_cur_gb, const  LSMLIB_REAL *rhs, const  int *ilo_rhs_gb, const  int *ihi_rhs_gb, const  int *jlo_rhs_gb, const  int *jhi_rhs_gb, const  int *klo_rhs_gb, const  int *khi_rhs_gb, const  LSMLIB_REAL *dt, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_utilities1d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/utilities/</path>
    <filename>lsm__utilities1d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_MAX_NORM_DIFF</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_COMPUTE_STABLE_ADVECTION_DT</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_COMPUTE_STABLE_NORMAL_VEL_DT</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_SURFACE_INTEGRAL</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_MAX_NORM_DIFF_CONTROL_VOLUME</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_COMPUTE_STABLE_ADVECTION_DT_CONTROL_VOLUME</name>
      <anchor>a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME</name>
      <anchor>a8</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a9</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a10</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM1D_SURFACE_INTEGRAL_CONTROL_VOLUME</name>
      <anchor>a11</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_MAX_NORM_DIFF</name>
      <anchor>a12</anchor>
      <arglist>(LSMLIB_REAL *max_norm_diff, const  LSMLIB_REAL *field1, const  int *ilo_field1_gb, const  int *ihi_field1_gb, const  LSMLIB_REAL *field2, const  int *ilo_field2_gb, const  int *ihi_field2_gb, const  int *ilo_ib, const  int *ihi_ib)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_COMPUTE_STABLE_ADVECTION_DT</name>
      <anchor>a13</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_x, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *ilo_ib, const  int *ihi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *cfl_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_COMPUTE_STABLE_NORMAL_VEL_DT</name>
      <anchor>a14</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_n, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  LSMLIB_REAL *phi_x_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *ilo_ib, const  int *ihi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *cfl_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO</name>
      <anchor>a15</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO</name>
      <anchor>a16</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_SURFACE_INTEGRAL</name>
      <anchor>a17</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  LSMLIB_REAL *phi_x, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_MAX_NORM_DIFF_CONTROL_VOLUME</name>
      <anchor>a18</anchor>
      <arglist>(LSMLIB_REAL *max_norm_diff, const  LSMLIB_REAL *field1, const  int *ilo_field1_gb, const  int *ihi_field1_gb, const  LSMLIB_REAL *field2, const  int *ilo_field2_gb, const  int *ihi_field2_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_COMPUTE_STABLE_ADVECTION_DT_CONTROL_VOLUME</name>
      <anchor>a19</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_x, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *cfl_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME</name>
      <anchor>a20</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_n, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  LSMLIB_REAL *phi_x_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *cfl_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a21</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a22</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM1D_SURFACE_INTEGRAL_CONTROL_VOLUME</name>
      <anchor>a23</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  LSMLIB_REAL *phi_x, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_utilities2d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/utilities/</path>
    <filename>lsm__utilities2d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_MAX_NORM_DIFF</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_AVE_ABS_DIFF</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_VOXEL_COUNT_GREATER_THAN_ZERO</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_VOXEL_COUNT_LESS_THAN_ZERO</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_STABLE_ADVECTION_DT</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_STABLE_NORMAL_VEL_DT</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO</name>
      <anchor>a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO</name>
      <anchor>a8</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_SURFACE_INTEGRAL</name>
      <anchor>a9</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_SURFACE_INTEGRAL_DELTA</name>
      <anchor>a10</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_MAX_NORM_DIFF_CONTROL_VOLUME</name>
      <anchor>a11</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_STABLE_ADVECTION_DT_CONTROL_VOLUME</name>
      <anchor>a12</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME</name>
      <anchor>a13</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT_CONTROL_VOLUME</name>
      <anchor>a14</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a15</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a16</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_SURFACE_INTEGRAL_CONTROL_VOLUME</name>
      <anchor>a17</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_VOXEL_COUNT_GREATER_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a18</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_VOXEL_COUNT_LESS_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a19</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_SURFACE_INTEGRAL_DELTA_CONTROL_VOLUME</name>
      <anchor>a20</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_MAX_NORM_DIFF</name>
      <anchor>a21</anchor>
      <arglist>(LSMLIB_REAL *max_norm_diff, const  LSMLIB_REAL *field1, const  int *ilo_field1_gb, const  int *ihi_field1_gb, const  int *jlo_field1_gb, const  int *jhi_field1_gb, const  LSMLIB_REAL *field2, const  int *ilo_field2_gb, const  int *ihi_field2_gb, const  int *jlo_field2_gb, const  int *jhi_field2_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_AVE_ABS_DIFF</name>
      <anchor>a22</anchor>
      <arglist>(LSMLIB_REAL *ave_abs_diff, const  LSMLIB_REAL *field1, const  int *ilo_field1_gb, const  int *ihi_field1_gb, const  int *jlo_field1_gb, const  int *jhi_field1_gb, const  LSMLIB_REAL *field2, const  int *ilo_field2_gb, const  int *ihi_field2_gb, const  int *jlo_field2_gb, const  int *jhi_field2_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_VOXEL_COUNT_GREATER_THAN_ZERO</name>
      <anchor>a23</anchor>
      <arglist>(int *count, const  LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_VOXEL_COUNT_LESS_THAN_ZERO</name>
      <anchor>a24</anchor>
      <arglist>(int *count, const  LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_STABLE_ADVECTION_DT</name>
      <anchor>a25</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *cfl_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_STABLE_NORMAL_VEL_DT</name>
      <anchor>a26</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_n, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *cfl_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT</name>
      <anchor>a27</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_n, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *cfl_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO</name>
      <anchor>a28</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  int *jlo_F_gb, const  int *jhi_F_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO</name>
      <anchor>a29</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  int *jlo_F_gb, const  int *jhi_F_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_SURFACE_INTEGRAL</name>
      <anchor>a30</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  int *jlo_F_gb, const  int *jhi_F_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_SURFACE_INTEGRAL_DELTA</name>
      <anchor>a31</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  int *jlo_F_gb, const  int *jhi_F_gb, const  LSMLIB_REAL *delta_phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_MAX_NORM_DIFF_CONTROL_VOLUME</name>
      <anchor>a32</anchor>
      <arglist>(LSMLIB_REAL *max_norm_diff, const  LSMLIB_REAL *field1, const  int *ilo_field1_gb, const  int *ihi_field1_gb, const  int *jlo_field1_gb, const  int *jhi_field1_gb, const  LSMLIB_REAL *field2, const  int *ilo_field2_gb, const  int *ihi_field2_gb, const  int *jlo_field2_gb, const  int *jhi_field2_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_STABLE_ADVECTION_DT_CONTROL_VOLUME</name>
      <anchor>a33</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *cfl_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME</name>
      <anchor>a34</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_n, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *cfl_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT_CONTROL_VOLUME</name>
      <anchor>a35</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_n, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *cfl_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a36</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  int *jlo_F_gb, const  int *jhi_F_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a37</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  int *jlo_F_gb, const  int *jhi_F_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_SURFACE_INTEGRAL_CONTROL_VOLUME</name>
      <anchor>a38</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  int *jlo_F_gb, const  int *jhi_F_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_VOXEL_COUNT_GREATER_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a39</anchor>
      <arglist>(int *count, const  LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_VOXEL_COUNT_LESS_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a40</anchor>
      <arglist>(int *count, const  LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_SURFACE_INTEGRAL_DELTA_CONTROL_VOLUME</name>
      <anchor>a41</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  int *jlo_F_gb, const  int *jhi_F_gb, const  LSMLIB_REAL *delta_phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  LSMLIB_REAL *grad_phi_mag, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_utilities2d_local.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/utilities/</path>
    <filename>lsm__utilities2d__local_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_MAX_NORM_DIFF_LOCAL</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_AVE_ABS_DIFF_LOCAL</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_STABLE_ADVECTION_DT_LOCAL</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_STABLE_NORMAL_VEL_DT_LOCAL</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM2D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT_LOCAL</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_MAX_NORM_DIFF_LOCAL</name>
      <anchor>a5</anchor>
      <arglist>(LSMLIB_REAL *max_norm_diff, const  LSMLIB_REAL *field1, const  int *ilo_field1_gb, const  int *ihi_field1_gb, const  int *jlo_field1_gb, const  int *jhi_field1_gb, const  LSMLIB_REAL *field2, const  int *ilo_field2_gb, const  int *ihi_field2_gb, const  int *jlo_field2_gb, const  int *jhi_field2_gb, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_AVE_ABS_DIFF_LOCAL</name>
      <anchor>a6</anchor>
      <arglist>(LSMLIB_REAL *ave_abs_diff, const  LSMLIB_REAL *field1, const  int *ilo_field1_gb, const  int *ihi_field1_gb, const  int *jlo_field1_gb, const  int *jhi_field1_gb, const  LSMLIB_REAL *field2, const  int *ilo_field2_gb, const  int *ihi_field2_gb, const  int *jlo_field2_gb, const  int *jhi_field2_gb, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_STABLE_ADVECTION_DT_LOCAL</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *cfl_number, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_STABLE_NORMAL_VEL_DT_LOCAL</name>
      <anchor>a8</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_n, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *cfl_number, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM2D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT_LOCAL</name>
      <anchor>a9</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_n, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *cfl_number, const  int *index_x, const  int *index_y, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_utilities3d.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/utilities/</path>
    <filename>lsm__utilities3d_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_MAX_NORM_DIFF</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_AVE_ABS_DIFF</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_VOXEL_COUNT_GREATER_THAN_ZERO</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_VOXEL_COUNT_LESS_THAN_ZERO</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_STABLE_ADVECTION_DT</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_STABLE_NORMAL_VEL_DT</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO</name>
      <anchor>a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO</name>
      <anchor>a8</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_SURFACE_INTEGRAL</name>
      <anchor>a9</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_MAX_NORM_DIFF_CONTROL_VOLUME</name>
      <anchor>a10</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_STABLE_ADVECTION_DT_CONTROL_VOLUME</name>
      <anchor>a11</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME</name>
      <anchor>a12</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT_CONTROL_VOLUME</name>
      <anchor>a13</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a14</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a15</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_SURFACE_INTEGRAL_CONTROL_VOLUME</name>
      <anchor>a16</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_VOXEL_COUNT_GREATER_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a17</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_VOXEL_COUNT_LESS_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a18</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_MAX_NORM_DIFF</name>
      <anchor>a19</anchor>
      <arglist>(LSMLIB_REAL *max_norm_diff, const  LSMLIB_REAL *field1, const  int *ilo_field1_gb, const  int *ihi_field1_gb, const  int *jlo_field1_gb, const  int *jhi_field1_gb, const  int *klo_field1_gb, const  int *khi_field1_gb, const  LSMLIB_REAL *field2, const  int *ilo_field2_gb, const  int *ihi_field2_gb, const  int *jlo_field2_gb, const  int *jhi_field2_gb, const  int *klo_field2_gb, const  int *khi_field2_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_AVE_ABS_DIFF</name>
      <anchor>a20</anchor>
      <arglist>(LSMLIB_REAL *ave_abs_diff, const  LSMLIB_REAL *field1, const  int *ilo_field1_gb, const  int *ihi_field1_gb, const  int *jlo_field1_gb, const  int *jhi_field1_gb, const  int *klo_field1_gb, const  int *khi_field1_gb, const  LSMLIB_REAL *field2, const  int *ilo_field2_gb, const  int *ihi_field2_gb, const  int *jlo_field2_gb, const  int *jhi_field2_gb, const  int *klo_field2_gb, const  int *khi_field2_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_VOXEL_COUNT_GREATER_THAN_ZERO</name>
      <anchor>a21</anchor>
      <arglist>(int *count, const  LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *klo_gb, const  int *khi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_VOXEL_COUNT_LESS_THAN_ZERO</name>
      <anchor>a22</anchor>
      <arglist>(int *count, const  LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *klo_gb, const  int *khi_gb, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_STABLE_ADVECTION_DT</name>
      <anchor>a23</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  LSMLIB_REAL *vel_z, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *klo_vel_gb, const  int *khi_vel_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *cfl_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_STABLE_NORMAL_VEL_DT</name>
      <anchor>a24</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_n, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *klo_vel_gb, const  int *khi_vel_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *cfl_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT</name>
      <anchor>a25</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_n, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *cfl_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO</name>
      <anchor>a26</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  int *jlo_F_gb, const  int *jhi_F_gb, const  int *klo_F_gb, const  int *khi_F_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO</name>
      <anchor>a27</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  int *jlo_F_gb, const  int *jhi_F_gb, const  int *klo_F_gb, const  int *khi_F_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_SURFACE_INTEGRAL</name>
      <anchor>a28</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  int *jlo_F_gb, const  int *jhi_F_gb, const  int *klo_F_gb, const  int *khi_F_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_MAX_NORM_DIFF_CONTROL_VOLUME</name>
      <anchor>a29</anchor>
      <arglist>(LSMLIB_REAL *max_norm_diff, const  LSMLIB_REAL *field1, const  int *ilo_field1_gb, const  int *ihi_field1_gb, const  int *jlo_field1_gb, const  int *jhi_field1_gb, const  int *klo_field1_gb, const  int *khi_field1_gb, const  LSMLIB_REAL *field2, const  int *ilo_field2_gb, const  int *ihi_field2_gb, const  int *jlo_field2_gb, const  int *jhi_field2_gb, const  int *klo_field2_gb, const  int *khi_field2_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *klo_control_vol_gb, const  int *khi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_STABLE_ADVECTION_DT_CONTROL_VOLUME</name>
      <anchor>a30</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  LSMLIB_REAL *vel_z, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *klo_vel_gb, const  int *khi_vel_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *klo_control_vol_gb, const  int *khi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *cfl_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME</name>
      <anchor>a31</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_n, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *klo_vel_gb, const  int *khi_vel_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *klo_control_vol_gb, const  int *khi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *cfl_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT_CONTROL_VOLUME</name>
      <anchor>a32</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_n, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *klo_control_vol_gb, const  int *khi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *cfl_number)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_VOLUME_INTEGRAL_PHI_LESS_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a33</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  int *jlo_F_gb, const  int *jhi_F_gb, const  int *klo_F_gb, const  int *khi_F_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *klo_control_vol_gb, const  int *khi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_VOLUME_INTEGRAL_PHI_GREATER_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a34</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  int *jlo_F_gb, const  int *jhi_F_gb, const  int *klo_F_gb, const  int *khi_F_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *klo_control_vol_gb, const  int *khi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_SURFACE_INTEGRAL_CONTROL_VOLUME</name>
      <anchor>a35</anchor>
      <arglist>(LSMLIB_REAL *int_F, const  LSMLIB_REAL *F, const  int *ilo_F_gb, const  int *ihi_F_gb, const  int *jlo_F_gb, const  int *jhi_F_gb, const  int *klo_F_gb, const  int *khi_F_gb, const  LSMLIB_REAL *phi, const  int *ilo_phi_gb, const  int *ihi_phi_gb, const  int *jlo_phi_gb, const  int *jhi_phi_gb, const  int *klo_phi_gb, const  int *khi_phi_gb, const  LSMLIB_REAL *phi_x, const  LSMLIB_REAL *phi_y, const  LSMLIB_REAL *phi_z, const  int *ilo_grad_phi_gb, const  int *ihi_grad_phi_gb, const  int *jlo_grad_phi_gb, const  int *jhi_grad_phi_gb, const  int *klo_grad_phi_gb, const  int *khi_grad_phi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *klo_control_vol_gb, const  int *khi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_ib, const  int *ihi_ib, const  int *jlo_ib, const  int *jhi_ib, const  int *klo_ib, const  int *khi_ib, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *epsilon)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_VOXEL_COUNT_GREATER_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a36</anchor>
      <arglist>(int *count, const  LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *klo_gb, const  int *khi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *klo_control_vol_gb, const  int *khi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_VOXEL_COUNT_LESS_THAN_ZERO_CONTROL_VOLUME</name>
      <anchor>a37</anchor>
      <arglist>(int *count, const  LSMLIB_REAL *phi, const  int *ilo_gb, const  int *ihi_gb, const  int *jlo_gb, const  int *jhi_gb, const  int *klo_gb, const  int *khi_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *klo_control_vol_gb, const  int *khi_control_vol_gb, const  int *control_vol_sgn, const  int *ilo_fb, const  int *ihi_fb, const  int *jlo_fb, const  int *jhi_fb, const  int *klo_fb, const  int *khi_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsm_utilities3d_local.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/utilities/</path>
    <filename>lsm__utilities3d__local_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_MAX_NORM_DIFF_LOCAL</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_STABLE_ADVECTION_DT_LOCAL</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_STABLE_NORMAL_VEL_DT_LOCAL</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME_LOCAL</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM3D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT_LOCAL</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_MAX_NORM_DIFF_LOCAL</name>
      <anchor>a5</anchor>
      <arglist>(LSMLIB_REAL *max_norm_diff, const  LSMLIB_REAL *field1, const  int *ilo_field1_gb, const  int *ihi_field1_gb, const  int *jlo_field1_gb, const  int *jhi_field1_gb, const  int *klo_field1_gb, const  int *khi_field1_gb, const  LSMLIB_REAL *field2, const  int *ilo_field2_gb, const  int *ihi_field2_gb, const  int *jlo_field2_gb, const  int *jhi_field2_gb, const  int *klo_field2_gb, const  int *khi_field2_gb, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_STABLE_ADVECTION_DT_LOCAL</name>
      <anchor>a6</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_x, const  LSMLIB_REAL *vel_y, const  LSMLIB_REAL *vel_z, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *klo_vel_gb, const  int *khi_vel_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *cfl_number, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_STABLE_NORMAL_VEL_DT_LOCAL</name>
      <anchor>a7</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_n, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *klo_vel_gb, const  int *khi_vel_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *cfl_number, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_STABLE_NORMAL_VEL_DT_CONTROL_VOLUME_LOCAL</name>
      <anchor>a8</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_n, const  int *ilo_vel_gb, const  int *ihi_vel_gb, const  int *jlo_vel_gb, const  int *jhi_vel_gb, const  int *klo_vel_gb, const  int *khi_vel_gb, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *control_vol, const  int *ilo_control_vol_gb, const  int *ihi_control_vol_gb, const  int *jlo_control_vol_gb, const  int *jhi_control_vol_gb, const  int *klo_control_vol_gb, const  int *khi_control_vol_gb, const  int *control_vol_sgn, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *cfl_number, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>LSM3D_COMPUTE_STABLE_CONST_NORMAL_VEL_DT_LOCAL</name>
      <anchor>a9</anchor>
      <arglist>(LSMLIB_REAL *dt, const  LSMLIB_REAL *vel_n, const  LSMLIB_REAL *phi_x_plus, const  LSMLIB_REAL *phi_y_plus, const  LSMLIB_REAL *phi_z_plus, const  int *ilo_grad_phi_plus_gb, const  int *ihi_grad_phi_plus_gb, const  int *jlo_grad_phi_plus_gb, const  int *jhi_grad_phi_plus_gb, const  int *klo_grad_phi_plus_gb, const  int *khi_grad_phi_plus_gb, const  LSMLIB_REAL *phi_x_minus, const  LSMLIB_REAL *phi_y_minus, const  LSMLIB_REAL *phi_z_minus, const  int *ilo_grad_phi_minus_gb, const  int *ihi_grad_phi_minus_gb, const  int *jlo_grad_phi_minus_gb, const  int *jhi_grad_phi_minus_gb, const  int *klo_grad_phi_minus_gb, const  int *khi_grad_phi_minus_gb, const  LSMLIB_REAL *dx, const  LSMLIB_REAL *dy, const  LSMLIB_REAL *dz, const  LSMLIB_REAL *cfl_number, const  int *index_x, const  int *index_y, const  int *index_z, const  int *nlo_index, const  int *nhi_index, const  unsigned char *narrow_band, const  int *ilo_nb_gb, const  int *ihi_nb_gb, const  int *jlo_nb_gb, const  int *jhi_nb_gb, const  int *klo_nb_gb, const  int *khi_nb_gb, const  unsigned char *mark_fb)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>LSMLIB_DefaultParameters.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/parallel/</path>
    <filename>LSMLIB__DefaultParameters_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>LSM_DIM_MAX</name>
      <anchor>a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM_DEFAULT_SPATIAL_DERIVATIVE_TYPE</name>
      <anchor>a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM_DEFAULT_SPATIAL_DERIVATIVE_WENO_ORDER</name>
      <anchor>a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM_DEFAULT_SPATIAL_DERIVATIVE_ENO_ORDER</name>
      <anchor>a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM_DEFAULT_TVD_RUNGE_KUTTA_ORDER</name>
      <anchor>a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM_DEFAULT_CFL_NUMBER</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>LSM_DEFAULT_VERBOSE_MODE</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>lsmlib_matlab_toolbox.dox</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/doc/doxygen/users_guide/</path>
    <filename>lsmlib__matlab__toolbox_8dox</filename>
  </compound>
  <compound kind="file">
    <name>lsmlib_toolbox_package.dox</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/doc/doxygen/users_guide/</path>
    <filename>lsmlib__toolbox__package_8dox</filename>
  </compound>
  <compound kind="file">
    <name>mainpage.dox</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/</path>
    <filename>mainpage_8dox</filename>
  </compound>
  <compound kind="file">
    <name>toolbox/manual.dox</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/toolbox/</path>
    <filename>toolbox_2manual_8dox</filename>
    <namespace>LSMLIB</namespace>
  </compound>
  <compound kind="file">
    <name>serial/manual.dox</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/serial/</path>
    <filename>serial_2manual_8dox</filename>
    <namespace>LSMLIB</namespace>
  </compound>
  <compound kind="file">
    <name>parallel/manual.dox</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/parallel/</path>
    <filename>parallel_2manual_8dox</filename>
    <namespace>LSMLIB</namespace>
  </compound>
  <compound kind="file">
    <name>matlab/manual.dox</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/matlab/</path>
    <filename>matlab_2manual_8dox</filename>
    <namespace>LSMLIB</namespace>
  </compound>
  <compound kind="file">
    <name>OrthogonalizationAlgorithm.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/parallel/</path>
    <filename>OrthogonalizationAlgorithm_8h</filename>
    <includes id="LevelSetMethodToolbox_8h" name="LevelSetMethodToolbox.h" local="yes">LevelSetMethodToolbox.h</includes>
    <includes id="FieldExtensionAlgorithm_8h" name="FieldExtensionAlgorithm.h" local="yes">FieldExtensionAlgorithm.h</includes>
    <includes id="BoundaryConditionModule_8h" name="BoundaryConditionModule.h" local="yes">BoundaryConditionModule.h</includes>
    <namespace>LSMLIB</namespace>
  </compound>
  <compound kind="file">
    <name>overview.dox</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/doc/doxygen/users_guide/</path>
    <filename>overview_8dox</filename>
  </compound>
  <compound kind="file">
    <name>parallel_lsmlib_package.dox</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/doc/doxygen/users_guide/</path>
    <filename>parallel__lsmlib__package_8dox</filename>
  </compound>
  <compound kind="file">
    <name>ReinitializationAlgorithm.h</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/src/parallel/</path>
    <filename>ReinitializationAlgorithm_8h</filename>
    <includes id="BoundaryConditionModule_8h" name="BoundaryConditionModule.h" local="yes">BoundaryConditionModule.h</includes>
    <includes id="LevelSetMethodToolbox_8h" name="LevelSetMethodToolbox.h" local="yes">LevelSetMethodToolbox.h</includes>
    <namespace>LSMLIB</namespace>
  </compound>
  <compound kind="file">
    <name>sample_input_file.dox</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/doc/doxygen/</path>
    <filename>sample__input__file_8dox</filename>
  </compound>
  <compound kind="file">
    <name>serial_lsmlib_package.dox</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/doc/doxygen/users_guide/</path>
    <filename>serial__lsmlib__package_8dox</filename>
  </compound>
  <compound kind="file">
    <name>users_guide.dox</name>
    <path>/Users/ktchu/Work/SerendipityResearch/Projects/scientific_software/release/lsmlib/doc/doxygen/users_guide/</path>
    <filename>users__guide_8dox</filename>
  </compound>
  <compound kind="struct">
    <name>_Grid</name>
    <filename>struct__Grid.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>num_dims</name>
      <anchor>o0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL</type>
      <name>x_lo</name>
      <anchor>o1</anchor>
      <arglist>[3]</arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL</type>
      <name>x_hi</name>
      <anchor>o2</anchor>
      <arglist>[3]</arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL</type>
      <name>x_lo_ghostbox</name>
      <anchor>o3</anchor>
      <arglist>[3]</arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL</type>
      <name>x_hi_ghostbox</name>
      <anchor>o4</anchor>
      <arglist>[3]</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>grid_dims</name>
      <anchor>o5</anchor>
      <arglist>[3]</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>grid_dims_ghostbox</name>
      <anchor>o6</anchor>
      <arglist>[3]</arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL</type>
      <name>dx</name>
      <anchor>o7</anchor>
      <arglist>[3]</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_gridpts</name>
      <anchor>o8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>ilo_gb</name>
      <anchor>o9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>ihi_gb</name>
      <anchor>o10</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>jlo_gb</name>
      <anchor>o11</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>jhi_gb</name>
      <anchor>o12</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>klo_gb</name>
      <anchor>o13</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>khi_gb</name>
      <anchor>o14</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>ilo_fb</name>
      <anchor>o15</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>ihi_fb</name>
      <anchor>o16</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>jlo_fb</name>
      <anchor>o17</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>jhi_fb</name>
      <anchor>o18</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>klo_fb</name>
      <anchor>o19</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>khi_fb</name>
      <anchor>o20</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>ilo_D1_fb</name>
      <anchor>o21</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>ihi_D1_fb</name>
      <anchor>o22</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>jlo_D1_fb</name>
      <anchor>o23</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>jhi_D1_fb</name>
      <anchor>o24</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>klo_D1_fb</name>
      <anchor>o25</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>khi_D1_fb</name>
      <anchor>o26</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>ilo_D2_fb</name>
      <anchor>o27</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>ihi_D2_fb</name>
      <anchor>o28</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>jlo_D2_fb</name>
      <anchor>o29</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>jhi_D2_fb</name>
      <anchor>o30</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>klo_D2_fb</name>
      <anchor>o31</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>khi_D2_fb</name>
      <anchor>o32</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>ilo_D3_fb</name>
      <anchor>o33</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>ihi_D3_fb</name>
      <anchor>o34</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>jlo_D3_fb</name>
      <anchor>o35</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>jhi_D3_fb</name>
      <anchor>o36</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>klo_D3_fb</name>
      <anchor>o37</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>khi_D3_fb</name>
      <anchor>o38</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_nb_levels</name>
      <anchor>o39</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned char</type>
      <name>mark_gb</name>
      <anchor>o40</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned char</type>
      <name>mark_D1</name>
      <anchor>o41</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned char</type>
      <name>mark_D2</name>
      <anchor>o42</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned char</type>
      <name>mark_D3</name>
      <anchor>o43</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned char</type>
      <name>mark_fb</name>
      <anchor>o44</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL</type>
      <name>beta</name>
      <anchor>o45</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL</type>
      <name>gamma</name>
      <anchor>o46</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>_LSM_DataArrays</name>
    <filename>struct__LSM__DataArrays.html</filename>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi</name>
      <anchor>o0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_stage1</name>
      <anchor>o1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_stage2</name>
      <anchor>o2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_next</name>
      <anchor>o3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi0</name>
      <anchor>o4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_prev</name>
      <anchor>o5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_extra</name>
      <anchor>o6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>mask</name>
      <anchor>o7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>lse_rhs</name>
      <anchor>o8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_x_plus</name>
      <anchor>o9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_x_minus</name>
      <anchor>o10</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_x</name>
      <anchor>o11</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_y_plus</name>
      <anchor>o12</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_y_minus</name>
      <anchor>o13</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_y</name>
      <anchor>o14</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_z_plus</name>
      <anchor>o15</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_z_minus</name>
      <anchor>o16</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_z</name>
      <anchor>o17</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>D1</name>
      <anchor>o18</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>D2</name>
      <anchor>o19</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>D3</name>
      <anchor>o20</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_xx</name>
      <anchor>o21</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_yy</name>
      <anchor>o22</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_xy</name>
      <anchor>o23</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_zz</name>
      <anchor>o24</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_xz</name>
      <anchor>o25</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>phi_yz</name>
      <anchor>o26</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>normal_velocity</name>
      <anchor>o27</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>external_velocity_x</name>
      <anchor>o28</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>external_velocity_y</name>
      <anchor>o29</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>external_velocity_z</name>
      <anchor>o30</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned char *</type>
      <name>narrow_band</name>
      <anchor>o31</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_index_pts</name>
      <anchor>o32</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>index_x</name>
      <anchor>o33</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>index_y</name>
      <anchor>o34</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>index_z</name>
      <anchor>o35</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>n_lo</name>
      <anchor>o36</anchor>
      <arglist>[10]</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>n_hi</name>
      <anchor>o37</anchor>
      <arglist>[10]</arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>index_outer_pts</name>
      <anchor>o38</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_alloc_index_outer_pts</name>
      <anchor>o39</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>nlo_outer_plus</name>
      <anchor>o40</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>nhi_outer_plus</name>
      <anchor>o41</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>nlo_outer_minus</name>
      <anchor>o42</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>nhi_outer_minus</name>
      <anchor>o43</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned char *</type>
      <name>solid_narrow_band</name>
      <anchor>o44</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>solid_num_index_pts</name>
      <anchor>o45</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>solid_index_x</name>
      <anchor>o46</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>solid_index_y</name>
      <anchor>o47</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>solid_index_z</name>
      <anchor>o48</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>solid_n_lo</name>
      <anchor>o49</anchor>
      <arglist>[10]</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>solid_n_hi</name>
      <anchor>o50</anchor>
      <arglist>[10]</arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>solid_normal_x</name>
      <anchor>o51</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>solid_normal_y</name>
      <anchor>o52</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL *</type>
      <name>solid_normal_z</name>
      <anchor>o53</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>HeapNode</name>
    <filename>structHeapNode.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>grid_idx</name>
      <anchor>o0</anchor>
      <arglist>[FMM_HEAP_MAX_NDIM]</arglist>
    </member>
    <member kind="variable">
      <type>LSMLIB_REAL</type>
      <name>value</name>
      <anchor>o1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>heap_pos</name>
      <anchor>o2</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>lsmlib_users_guide</name>
    <title>LSMLIB User&apos;s Guide</title>
    <filename>lsmlib_users_guide</filename>
  </compound>
  <compound kind="page">
    <name>users_guide_overview</name>
    <title>Overview</title>
    <filename>users_guide_overview</filename>
  </compound>
  <compound kind="page">
    <name>users_guide_lsmlib_toolbox_package</name>
    <title>LSMLIB Toolbox Package</title>
    <filename>users_guide_lsmlib_toolbox_package</filename>
    <docanchor>toolbox_notes</docanchor>
  </compound>
  <compound kind="page">
    <name>users_guide_serial_lsmlib_package</name>
    <title>Serial LSMLIB Package</title>
    <filename>users_guide_serial_lsmlib_package</filename>
    <docanchor>serial_boundary_conditions</docanchor>
    <docanchor>serial_fast_marching_method</docanchor>
    <docanchor>serial_initialization</docanchor>
    <docanchor>serial_grid_management</docanchor>
    <docanchor>serial_introduction</docanchor>
    <docanchor>serial_utilities</docanchor>
    <docanchor>serial_usage</docanchor>
  </compound>
  <compound kind="page">
    <name>users_guide_parallel_lsmlib_package</name>
    <title>Parallel LSMLIB Package</title>
    <filename>users_guide_parallel_lsmlib_package</filename>
    <docanchor>parallel_notes</docanchor>
    <docanchor>parallel_input</docanchor>
    <docanchor>parallel_velocity_field</docanchor>
    <docanchor>parallel_initialization</docanchor>
    <docanchor>parallel_introduction</docanchor>
    <docanchor>parallel_usage</docanchor>
    <docanchor>parallel_timestep_size</docanchor>
  </compound>
  <compound kind="page">
    <name>users_guide_lsmlib_matlab_toolbox</name>
    <title>LSMLIB MATLAB Toolbox</title>
    <filename>users_guide_lsmlib_matlab_toolbox</filename>
    <docanchor>matlab_fmm</docanchor>
    <docanchor>matlab_reinitialization</docanchor>
    <docanchor>matlab_introduction</docanchor>
    <docanchor>matlab_spatial_derivative</docanchor>
    <docanchor>matlab_tvdrk</docanchor>
    <docanchor>matlab_time_evolution</docanchor>
    <docanchor>matlab_usage</docanchor>
  </compound>
  <compound kind="page">
    <name>sample_input_file</name>
    <title>Sample Input File</title>
    <filename>sample_input_file</filename>
  </compound>
  <compound kind="page">
    <name>package_lsm_toolbox</name>
    <title>LSMLIB Toolbox Package</title>
    <filename>package_lsm_toolbox</filename>
  </compound>
  <compound kind="page">
    <name>package_lsm_serial</name>
    <title>Serial LSMLIB Package</title>
    <filename>package_lsm_serial</filename>
  </compound>
  <compound kind="page">
    <name>package_lsm_parallel</name>
    <title>Parallel LSMLIB Package</title>
    <filename>package_lsm_parallel</filename>
  </compound>
  <compound kind="page">
    <name>package_lsm_matlab</name>
    <title>LSMLIB MATLAB Toolbox</title>
    <filename>package_lsm_matlab</filename>
  </compound>
  <compound kind="namespace">
    <name>geom</name>
    <filename>namespacegeom.html</filename>
  </compound>
  <compound kind="namespace">
    <name>hier</name>
    <filename>namespacehier.html</filename>
  </compound>
  <compound kind="namespace">
    <name>LSMLIB</name>
    <filename>namespaceLSMLIB.html</filename>
    <class kind="class">LSMLIB::BoundaryConditionModule</class>
    <class kind="class">LSMLIB::FieldExtensionAlgorithm</class>
    <class kind="class">LSMLIB::LevelSetFunctionIntegrator</class>
    <class kind="class">LSMLIB::LevelSetFunctionIntegratorStrategy</class>
    <class kind="class">LSMLIB::LevelSetMethodAlgorithm</class>
    <class kind="class">LSMLIB::LevelSetMethodGriddingAlgorithm</class>
    <class kind="class">LSMLIB::LevelSetMethodGriddingStrategy</class>
    <class kind="class">LSMLIB::LevelSetMethodPatchStrategy</class>
    <class kind="class">LSMLIB::LevelSetMethodToolbox</class>
    <class kind="class">LSMLIB::LevelSetMethodVelocityFieldStrategy</class>
    <class kind="class">LSMLIB::OrthogonalizationAlgorithm</class>
    <class kind="class">LSMLIB::ReinitializationAlgorithm</class>
    <member kind="enumeration">
      <name>LEVEL_SET_FCN_TYPE</name>
      <anchor>a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>PHI</name>
      <anchor>a5a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>PSI</name>
      <anchor>a5a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumeration">
      <name>SPATIAL_DERIVATIVE_TYPE</name>
      <anchor>a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>ENO</name>
      <anchor>a6a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>WENO</name>
      <anchor>a6a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>UNKNOWN</name>
      <anchor>a6a4</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>LSMLIB::BoundaryConditionModule</name>
    <filename>classLSMLIB_1_1BoundaryConditionModule.html</filename>
    <templarg>DIM</templarg>
    <member kind="enumeration">
      <name>BOUNDARY_CONDITION_TYPE</name>
      <anchor>w5</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>NONE</name>
      <anchor>w5w0</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>HOMOGENEOUS_NEUMANN</name>
      <anchor>w5w1</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>LINEAR_EXTRAPOLATION</name>
      <anchor>w5w2</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>SIGNED_LINEAR_EXTRAPOLATION</name>
      <anchor>w5w3</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>ANTI_PERIODIC</name>
      <anchor>w5w4</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BoundaryConditionModule</name>
      <anchor>z3_0</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; patch_hierarchy, const  IntVector&lt; DIM &gt; &amp;ghostcell_width)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BoundaryConditionModule</name>
      <anchor>z3_1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BoundaryConditionModule</name>
      <anchor>z3_2</anchor>
      <arglist>(const  BoundaryConditionModule&lt; DIM &gt; &amp;rhs)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~BoundaryConditionModule</name>
      <anchor>z3_3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>imposeBoundaryConditions</name>
      <anchor>z5_0</anchor>
      <arglist>(const  int phi_handle, const  IntVector&lt; DIM &gt; &amp;lower_bc, const  IntVector&lt; DIM &gt; &amp;upper_bc, const  SPATIAL_DERIVATIVE_TYPE spatial_derivative_type, const  int spatial_derivative_order, const  int component=-1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>imposeBoundaryConditionsOnPatch</name>
      <anchor>z5_1</anchor>
      <arglist>(Patch&lt; DIM &gt; &amp;patch, const  int phi_handle, const  IntVector&lt; DIM &gt; &amp;lower_bc, const  IntVector&lt; DIM &gt; &amp;upper_bc, const  SPATIAL_DERIVATIVE_TYPE spatial_derivative_type, const  int spatial_derivative_order, const  int component=-1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>imposeAntiPeriodicBCs</name>
      <anchor>z5_2</anchor>
      <arglist>(const  int phi_handle, const  IntVector&lt; DIM &gt; &amp;lower_bc, const  IntVector&lt; DIM &gt; &amp;upper_bc, const  int component=-1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>imposeAntiPeriodicBCsOnPatch</name>
      <anchor>z5_3</anchor>
      <arglist>(Patch&lt; DIM &gt; &amp;patch, const  int phi_handle, const  IntVector&lt; DIM &gt; &amp;lower_bc, const  IntVector&lt; DIM &gt; &amp;upper_bc, const  int component=-1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>imposeHomogeneousNeumannBCs</name>
      <anchor>z5_4</anchor>
      <arglist>(const  int phi_handle, const  IntVector&lt; DIM &gt; &amp;lower_bc, const  IntVector&lt; DIM &gt; &amp;upper_bc, const  SPATIAL_DERIVATIVE_TYPE spatial_derivative_type, const  int spatial_derivative_order, const  int component=-1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>imposeHomogeneousNeumannBCsOnPatch</name>
      <anchor>z5_5</anchor>
      <arglist>(Patch&lt; DIM &gt; &amp;patch, const  int phi_handle, const  IntVector&lt; DIM &gt; &amp;lower_bc, const  IntVector&lt; DIM &gt; &amp;upper_bc, const  SPATIAL_DERIVATIVE_TYPE spatial_derivative_type, const  int spatial_derivative_order, const  int component=-1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>imposeLinearExtrapolationBCs</name>
      <anchor>z5_6</anchor>
      <arglist>(const  int phi_handle, const  IntVector&lt; DIM &gt; &amp;lower_bc, const  IntVector&lt; DIM &gt; &amp;upper_bc, const  int component=-1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>imposeLinearExtrapolationBCsOnPatch</name>
      <anchor>z5_7</anchor>
      <arglist>(Patch&lt; DIM &gt; &amp;patch, const  int phi_handle, const  IntVector&lt; DIM &gt; &amp;lower_bc, const  IntVector&lt; DIM &gt; &amp;upper_bc, const  int component=-1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>imposeSignedLinearExtrapolationBCs</name>
      <anchor>z5_8</anchor>
      <arglist>(const  int phi_handle, const  IntVector&lt; DIM &gt; &amp;lower_bc, const  IntVector&lt; DIM &gt; &amp;upper_bc, const  int component=-1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>imposeSignedLinearExtrapolationBCsOnPatch</name>
      <anchor>z5_9</anchor>
      <arglist>(Patch&lt; DIM &gt; &amp;patch, const  int phi_handle, const  IntVector&lt; DIM &gt; &amp;lower_bc, const  IntVector&lt; DIM &gt; &amp;upper_bc, const  int component=-1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>resetHierarchyConfiguration</name>
      <anchor>z7_0</anchor>
      <arglist>(const  Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; patch_hierarchy, const  int coarsest_level, const  int finest_level, const  IntVector&lt; DIM &gt; &amp;ghostcell_width)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const  BoundaryConditionModule &amp;</type>
      <name>operator=</name>
      <anchor>z9_0</anchor>
      <arglist>(const  BoundaryConditionModule&lt; DIM &gt; &amp;rhs)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>void</type>
      <name>computeIndexSpaceOfNearestGhostLayer</name>
      <anchor>z11_0</anchor>
      <arglist>(IntVector&lt; DIM &gt; &amp;nearest_ghost_layer_lo, IntVector&lt; DIM &gt; &amp;nearest_ghost_layer_hi, const  int bdry_type, const  int bdry_location_idx, const  Box&lt; DIM &gt; &amp;fillbox)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>int</type>
      <name>computeIndexOffset</name>
      <anchor>z11_1</anchor>
      <arglist>(const  int bdry_type, const  int bdry_location_idx, const  Box&lt; DIM &gt; &amp;ghostbox)</arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt;</type>
      <name>d_patch_hierarchy</name>
      <anchor>p0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>IntVector&lt; DIM &gt;</type>
      <name>d_ghostcell_width</name>
      <anchor>p1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>IntVector&lt; DIM &gt;</type>
      <name>d_geom_periodic_dirs</name>
      <anchor>p2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Array&lt; Array&lt; Array&lt; BoundaryBox&lt; DIM &gt; &gt; &gt; &gt;</type>
      <name>d_boundary_boxes</name>
      <anchor>p3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Array&lt; Array&lt; bool &gt; &gt;</type>
      <name>d_touches_boundary</name>
      <anchor>p4</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>LSMLIB::FieldExtensionAlgorithm</name>
    <filename>classLSMLIB_1_1FieldExtensionAlgorithm.html</filename>
    <templarg>DIM</templarg>
    <member kind="function">
      <type></type>
      <name>FieldExtensionAlgorithm</name>
      <anchor>z13_0</anchor>
      <arglist>(Pointer&lt; Database &gt; input_db, Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int field_handle, const  int phi_handle, const  int control_volume_handle, const  string &amp;object_name=&quot;FieldExtensionAlgorithm&quot;, const  IntVector&lt; DIM &gt; &amp;phi_ghostcell_width=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>FieldExtensionAlgorithm</name>
      <anchor>z13_1</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int field_handle, const  int phi_handle, const  int control_volume_handle, const  SPATIAL_DERIVATIVE_TYPE spatial_derivative_type=ENO, const  int spatial_derivative_order=1, const  int tvd_runge_kutta_order=1, const  LSMLIB_REAL cfl_number=0.5, const  LSMLIB_REAL stop_distance=0.0, const  int max_iterations=0, const  LSMLIB_REAL iteration_stop_tolerance=0.0, const  bool verbose_mode=false, const  string &amp;object_name=&quot;FieldExtensionAlgorithm&quot;, const  IntVector&lt; DIM &gt; &amp;phi_ghostcell_width=0)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~FieldExtensionAlgorithm</name>
      <anchor>z13_2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>computeExtensionField</name>
      <anchor>z15_0</anchor>
      <arglist>(const  int phi_component=0, const  int max_iterations=-1, const  IntVector&lt; DIM &gt; &amp;lower_bc_phi=IntVector&lt; DIM &gt;(-1), const  IntVector&lt; DIM &gt; &amp;upper_bc_phi=IntVector&lt; DIM &gt;(-1), const  IntVector&lt; DIM &gt; &amp;lower_bc_ext=IntVector&lt; DIM &gt;(-1), const  IntVector&lt; DIM &gt; &amp;upper_bc_ext=IntVector&lt; DIM &gt;(-1))</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>computeExtensionFieldForSingleComponent</name>
      <anchor>z15_1</anchor>
      <arglist>(const  int component=0, const  int phi_component=0, const  int max_iterations=-1, const  IntVector&lt; DIM &gt; &amp;lower_bc_phi=IntVector&lt; DIM &gt;(-1), const  IntVector&lt; DIM &gt; &amp;upper_bc_phi=IntVector&lt; DIM &gt;(-1), const  IntVector&lt; DIM &gt; &amp;lower_bc_ext=IntVector&lt; DIM &gt;(-1), const  IntVector&lt; DIM &gt; &amp;upper_bc_ext=IntVector&lt; DIM &gt;(-1))</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>resetHierarchyConfiguration</name>
      <anchor>z17_0</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int coarsest_level, const  int finest_level)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>advanceFieldExtensionEqnUsingTVDRK1</name>
      <anchor>z19_0</anchor>
      <arglist>(const  LSMLIB_REAL dt, const  int field_component, const  int phi_component, const  IntVector&lt; DIM &gt; &amp;lower_bc_ext, const  IntVector&lt; DIM &gt; &amp;upper_bc_ext)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>advanceFieldExtensionEqnUsingTVDRK2</name>
      <anchor>z19_1</anchor>
      <arglist>(const  LSMLIB_REAL dt, const  int field_component, const  int phi_component, const  IntVector&lt; DIM &gt; &amp;lower_bc_ext, const  IntVector&lt; DIM &gt; &amp;upper_bc_ext)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>advanceFieldExtensionEqnUsingTVDRK3</name>
      <anchor>z19_2</anchor>
      <arglist>(const  LSMLIB_REAL dt, const  int field_component, const  int phi_component, const  IntVector&lt; DIM &gt; &amp;lower_bc_ext, const  IntVector&lt; DIM &gt; &amp;upper_bc_ext)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>computeFieldExtensionEqnRHS</name>
      <anchor>z19_3</anchor>
      <arglist>(const  int extension_field_handle, const  int phi_component)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>initializeVariables</name>
      <anchor>z21_0</anchor>
      <arglist>(const  IntVector&lt; DIM &gt; &amp;phi_ghostcell_width)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>initializeCommunicationObjects</name>
      <anchor>z21_1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>getFromInput</name>
      <anchor>z21_2</anchor>
      <arglist>(Pointer&lt; Database &gt; db)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>checkParameters</name>
      <anchor>z21_3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="public">
      <type>string</type>
      <name>d_object_name</name>
      <anchor>p0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>SPATIAL_DERIVATIVE_TYPE</type>
      <name>d_spatial_derivative_type</name>
      <anchor>p1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_spatial_derivative_order</name>
      <anchor>p2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_tvd_runge_kutta_order</name>
      <anchor>p3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LSMLIB_REAL</type>
      <name>d_cfl_number</name>
      <anchor>p4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LSMLIB_REAL</type>
      <name>d_stop_distance</name>
      <anchor>p5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_max_iterations</name>
      <anchor>p6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LSMLIB_REAL</type>
      <name>d_iteration_stop_tol</name>
      <anchor>p7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_verbose_mode</name>
      <anchor>p8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt;</type>
      <name>d_patch_hierarchy</name>
      <anchor>p9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; CartesianGridGeometry&lt; DIM &gt; &gt;</type>
      <name>d_grid_geometry</name>
      <anchor>p10</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_extension_field_handle</name>
      <anchor>p11</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_phi_handle</name>
      <anchor>p12</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_control_volume_handle</name>
      <anchor>p13</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>vector&lt; int &gt;</type>
      <name>d_extension_field_scr_handles</name>
      <anchor>p14</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>IntVector&lt; DIM &gt;</type>
      <name>d_ext_field_scratch_ghostcell_width</name>
      <anchor>p15</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_phi_scr_handle</name>
      <anchor>p16</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>IntVector&lt; DIM &gt;</type>
      <name>d_phi_scratch_ghostcell_width</name>
      <anchor>p17</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_rhs_handle</name>
      <anchor>p18</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_normal_vector_handle</name>
      <anchor>p19</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_grad_field_handle</name>
      <anchor>p20</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_grad_phi_plus_handle</name>
      <anchor>p21</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_grad_phi_minus_handle</name>
      <anchor>p22</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_stop_distance</name>
      <anchor>p23</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_max_iterations</name>
      <anchor>p24</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_iteration_stop_tol</name>
      <anchor>p25</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_hierarchy_configuration_needs_reset</name>
      <anchor>p26</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_num_field_components</name>
      <anchor>p27</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>ComponentSelector</type>
      <name>d_scratch_data</name>
      <anchor>p28</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; BoundaryConditionModule&lt; DIM &gt; &gt;</type>
      <name>d_phi_bc_module</name>
      <anchor>p29</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; BoundaryConditionModule&lt; DIM &gt; &gt;</type>
      <name>d_ext_field_bc_module</name>
      <anchor>p30</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Array&lt; Pointer&lt; RefineAlgorithm&lt; DIM &gt; &gt; &gt;</type>
      <name>d_extension_field_fill_bdry_alg</name>
      <anchor>p31</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Array&lt; Array&lt; Pointer&lt; RefineSchedule&lt; DIM &gt; &gt; &gt; &gt;</type>
      <name>d_extension_field_fill_bdry_sched</name>
      <anchor>p32</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; RefineAlgorithm&lt; DIM &gt; &gt;</type>
      <name>d_phi_fill_bdry_alg</name>
      <anchor>p33</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Array&lt; Pointer&lt; RefineSchedule&lt; DIM &gt; &gt; &gt;</type>
      <name>d_phi_fill_bdry_sched</name>
      <anchor>p34</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>LSMLIB::LevelSetFunctionIntegrator</name>
    <filename>classLSMLIB_1_1LevelSetFunctionIntegrator.html</filename>
    <templarg>DIM</templarg>
    <base>LSMLIB::LevelSetFunctionIntegratorStrategy</base>
    <member kind="function">
      <type></type>
      <name>LevelSetFunctionIntegrator</name>
      <anchor>z23_0</anchor>
      <arglist>(Pointer&lt; Database &gt; input_db, Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; patch_hierarchy, LevelSetMethodPatchStrategy&lt; DIM &gt; *lsm_patch_strategy, LevelSetMethodVelocityFieldStrategy&lt; DIM &gt; *lsm_velocity_field_strategy, const  int num_level_set_fcn_components=1, const  int codimension=1, const  string &amp;object_name=&quot;LevelSetFunctionIntegrator&quot;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LevelSetFunctionIntegrator</name>
      <anchor>z23_1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getPhiPatchDataHandle</name>
      <anchor>z25_0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getPsiPatchDataHandle</name>
      <anchor>z25_1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getControlVolumePatchDataHandle</name>
      <anchor>z25_2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual LSMLIB_REAL</type>
      <name>getStartTime</name>
      <anchor>z27_0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual LSMLIB_REAL</type>
      <name>getEndTime</name>
      <anchor>z27_1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual LSMLIB_REAL</type>
      <name>getCurrentTime</name>
      <anchor>z27_2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>endTimeReached</name>
      <anchor>z27_3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>numIntegrationStepsTaken</name>
      <anchor>z27_4</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>printClassData</name>
      <anchor>z27_5</anchor>
      <arglist>(ostream &amp;os) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getSpatialDerivativeType</name>
      <anchor>z28_0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getSpatialDerivativeOrder</name>
      <anchor>z28_1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getTVDRungeKuttaOrder</name>
      <anchor>z28_2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual LSMLIB_REAL</type>
      <name>computeStableDt</name>
      <anchor>z30_0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>advanceLevelSetFunctions</name>
      <anchor>z30_1</anchor>
      <arglist>(const  LSMLIB_REAL dt)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getReinitializationInterval</name>
      <anchor>z32_0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setReinitializationInterval</name>
      <anchor>z32_1</anchor>
      <arglist>(const  int interval)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>reinitializeLevelSetFunctions</name>
      <anchor>z32_2</anchor>
      <arglist>(const  LEVEL_SET_FCN_TYPE level_set_fcn=LSMLIB::PHI, const  int max_iterations=-1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getOrthogonalizationInterval</name>
      <anchor>z32_3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setOrthogonalizationInterval</name>
      <anchor>z32_4</anchor>
      <arglist>(const  int interval)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>orthogonalizeLevelSetFunctions</name>
      <anchor>z32_5</anchor>
      <arglist>(const  LEVEL_SET_FCN_TYPE level_set_fcn, const  int max_reinit_iterations=-1, const  int max_ortho_iterations=-1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>preprocessInitializeVelocityField</name>
      <anchor>z34_0</anchor>
      <arglist>(int &amp;phi_handle, int &amp;psi_handle, const  Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int level_number)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>postprocessInitializeVelocityField</name>
      <anchor>z34_1</anchor>
      <arglist>(const  Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int level_number)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>putToDatabase</name>
      <anchor>z36_0</anchor>
      <arglist>(Pointer&lt; Database &gt; db)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>initializeLevelData</name>
      <anchor>z38_0</anchor>
      <arglist>(const  Pointer&lt; BasePatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int level_number, const  double init_data_time, const  bool can_be_refined, const  bool initial_time, const  Pointer&lt; BasePatchLevel&lt; DIM &gt; &gt; old_level=Pointer&lt; BasePatchLevel&lt; DIM &gt; &gt;((0)), const  bool allocate_data=true)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setPhysicalBoundaryConditions</name>
      <anchor>z38_1</anchor>
      <arglist>(Patch&lt; DIM &gt; &amp;patch, const  double fill_time, const  IntVector&lt; DIM &gt; &amp;ghost_width_to_fill)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setBoundaryConditions</name>
      <anchor>z38_2</anchor>
      <arglist>(const  IntVector&lt; DIM &gt; &amp;lower_bc, const  IntVector&lt; DIM &gt; &amp;upper_bc, const  LEVEL_SET_FCN_TYPE level_set_fcn, const  int component=-1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>applyGradientDetector</name>
      <anchor>z40_0</anchor>
      <arglist>(const  Pointer&lt; BasePatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int level_number, const  double error_data_time, const  int tag_index, const  bool initial_time, const  bool uses_richardson_extrapolation_too)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>resetHierarchyConfiguration</name>
      <anchor>z40_1</anchor>
      <arglist>(Pointer&lt; BasePatchHierarchy&lt; DIM &gt; &gt; hierarchy, int coarsest_level, int finest_level)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual IntVector&lt; DIM &gt;</type>
      <name>getRefineOpStencilWidth</name>
      <anchor>z41_0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>preprocessRefine</name>
      <anchor>z41_1</anchor>
      <arglist>(Patch&lt; DIM &gt; &amp;fine, const  Patch&lt; DIM &gt; &amp;coarse, const  Box&lt; DIM &gt; &amp;fine_box, const  IntVector&lt; DIM &gt; &amp;ratio)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>postprocessRefine</name>
      <anchor>z41_2</anchor>
      <arglist>(Patch&lt; DIM &gt; &amp;fine, const  Patch&lt; DIM &gt; &amp;coarse, const  Box&lt; DIM &gt; &amp;fine_box, const  IntVector&lt; DIM &gt; &amp;ratio)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual IntVector&lt; DIM &gt;</type>
      <name>getCoarsenOpStencilWidth</name>
      <anchor>z41_3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>preprocessCoarsen</name>
      <anchor>z41_4</anchor>
      <arglist>(Patch&lt; DIM &gt; &amp;coarse, const  Patch&lt; DIM &gt; &amp;fine, const  Box&lt; DIM &gt; &amp;coarse_box, const  IntVector&lt; DIM &gt; &amp;ratio)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>postprocessCoarsen</name>
      <anchor>z41_5</anchor>
      <arglist>(Patch&lt; DIM &gt; &amp;coarse, const  Patch&lt; DIM &gt; &amp;fine, const  Box&lt; DIM &gt; &amp;coarse_box, const  IntVector&lt; DIM &gt; &amp;ratio)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>advanceLevelSetEqnUsingTVDRK1</name>
      <anchor>z43_0</anchor>
      <arglist>(const  LSMLIB_REAL dt)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>advanceLevelSetEqnUsingTVDRK2</name>
      <anchor>z43_1</anchor>
      <arglist>(const  LSMLIB_REAL dt)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>advanceLevelSetEqnUsingTVDRK3</name>
      <anchor>z43_2</anchor>
      <arglist>(const  LSMLIB_REAL dt)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>computeLevelSetEquationRHS</name>
      <anchor>z43_3</anchor>
      <arglist>(const  LEVEL_SET_FCN_TYPE level_set_fcn, const  int phi_handle, const  int component=0)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>addAdvectionTermToLevelSetEquationRHS</name>
      <anchor>z43_4</anchor>
      <arglist>(const  LEVEL_SET_FCN_TYPE level_set_fcn, const  int phi_handle, const  int component=0)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>addNormalVelocityTermToLevelSetEquationRHS</name>
      <anchor>z43_5</anchor>
      <arglist>(const  LEVEL_SET_FCN_TYPE level_set_fcn, const  int phi_handle, const  int component=0)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>initializeVariables</name>
      <anchor>z45_0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>initializeCommunicationObjects</name>
      <anchor>z45_1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>getFromInput</name>
      <anchor>z45_2</anchor>
      <arglist>(Pointer&lt; Database &gt; db, bool is_from_restart)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>getFromRestart</name>
      <anchor>z45_3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="public">
      <type>string</type>
      <name>d_object_name</name>
      <anchor>p0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_num_level_set_fcn_components</name>
      <anchor>p1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_codimension</name>
      <anchor>p2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LSMLIB_REAL</type>
      <name>d_start_time</name>
      <anchor>p3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LSMLIB_REAL</type>
      <name>d_end_time</name>
      <anchor>p4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LSMLIB_REAL</type>
      <name>d_cfl_number</name>
      <anchor>p5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>SPATIAL_DERIVATIVE_TYPE</type>
      <name>d_spatial_derivative_type</name>
      <anchor>p6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_spatial_derivative_order</name>
      <anchor>p7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_tvd_runge_kutta_order</name>
      <anchor>p8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_reinitialization_interval</name>
      <anchor>p9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LSMLIB_REAL</type>
      <name>d_reinitialization_stop_tol</name>
      <anchor>p10</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LSMLIB_REAL</type>
      <name>d_reinitialization_stop_dist</name>
      <anchor>p11</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_reinitialization_max_iters</name>
      <anchor>p12</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_orthogonalization_interval</name>
      <anchor>p13</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LSMLIB_REAL</type>
      <name>d_orthogonalization_stop_tol</name>
      <anchor>p14</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LSMLIB_REAL</type>
      <name>d_orthogonalization_stop_dist</name>
      <anchor>p15</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_orthogonalization_max_iters</name>
      <anchor>p16</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_AMR</name>
      <anchor>p17</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_regrid_interval</name>
      <anchor>p18</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_tag_buffer_width</name>
      <anchor>p19</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LSMLIB_REAL</type>
      <name>d_refinement_cutoff_value</name>
      <anchor>p20</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_verbose_mode</name>
      <anchor>p21</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LevelSetMethodPatchStrategy&lt; DIM &gt; *</type>
      <name>d_lsm_patch_strategy</name>
      <anchor>p22</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LevelSetMethodVelocityFieldStrategy&lt; DIM &gt; *</type>
      <name>d_lsm_velocity_field_strategy</name>
      <anchor>p23</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt;</type>
      <name>d_patch_hierarchy</name>
      <anchor>p24</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; CartesianGridGeometry&lt; DIM &gt; &gt;</type>
      <name>d_grid_geometry</name>
      <anchor>p25</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; ReinitializationAlgorithm&lt; DIM &gt; &gt;</type>
      <name>d_phi_reinitialization_alg</name>
      <anchor>p26</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; ReinitializationAlgorithm&lt; DIM &gt; &gt;</type>
      <name>d_psi_reinitialization_alg</name>
      <anchor>p27</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; OrthogonalizationAlgorithm&lt; DIM &gt; &gt;</type>
      <name>d_orthogonalization_alg</name>
      <anchor>p28</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>vector&lt; int &gt;</type>
      <name>d_phi_handles</name>
      <anchor>p29</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>vector&lt; int &gt;</type>
      <name>d_psi_handles</name>
      <anchor>p30</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_grad_phi_plus_handle</name>
      <anchor>p31</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_grad_psi_plus_handle</name>
      <anchor>p32</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_grad_phi_minus_handle</name>
      <anchor>p33</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_grad_psi_minus_handle</name>
      <anchor>p34</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_grad_phi_upwind_handle</name>
      <anchor>p35</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_grad_psi_upwind_handle</name>
      <anchor>p36</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_rhs_phi_handle</name>
      <anchor>p37</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_rhs_psi_handle</name>
      <anchor>p38</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_control_volume_handle</name>
      <anchor>p39</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>IntVector&lt; DIM &gt;</type>
      <name>d_level_set_ghostcell_width</name>
      <anchor>p40</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>ComponentSelector</type>
      <name>d_solution_variables</name>
      <anchor>p41</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>ComponentSelector</type>
      <name>d_time_advance_scratch_variables</name>
      <anchor>p42</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>ComponentSelector</type>
      <name>d_compute_stable_dt_scratch_variables</name>
      <anchor>p43</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>ComponentSelector</type>
      <name>d_reinitialization_scratch_variables</name>
      <anchor>p44</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>ComponentSelector</type>
      <name>d_orthogonalization_scratch_variables</name>
      <anchor>p45</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>ComponentSelector</type>
      <name>d_persistent_variables</name>
      <anchor>p46</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_reinitialization</name>
      <anchor>p47</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_reinitialization_stop_tol</name>
      <anchor>p48</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_reinitialization_stop_dist</name>
      <anchor>p49</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_reinitialization_max_iters</name>
      <anchor>p50</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_orthogonalization</name>
      <anchor>p51</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_orthogonalization_stop_tol</name>
      <anchor>p52</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_orthogonalization_stop_dist</name>
      <anchor>p53</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_orthogonalization_max_iters</name>
      <anchor>p54</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LSMLIB_REAL</type>
      <name>d_current_time</name>
      <anchor>p55</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_num_integration_steps_taken</name>
      <anchor>p56</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_reinitialization_count</name>
      <anchor>p57</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_orthogonalization_count</name>
      <anchor>p58</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LEVEL_SET_FCN_TYPE</type>
      <name>d_orthogonalization_evolved_field</name>
      <anchor>p59</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_regrid_count</name>
      <anchor>p60</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; BoundaryConditionModule&lt; DIM &gt; &gt;</type>
      <name>d_bc_module</name>
      <anchor>p61</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Array&lt; IntVector&lt; DIM &gt; &gt;</type>
      <name>d_lower_bc_phi</name>
      <anchor>p62</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Array&lt; IntVector&lt; DIM &gt; &gt;</type>
      <name>d_upper_bc_phi</name>
      <anchor>p63</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Array&lt; IntVector&lt; DIM &gt; &gt;</type>
      <name>d_lower_bc_psi</name>
      <anchor>p64</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Array&lt; IntVector&lt; DIM &gt; &gt;</type>
      <name>d_upper_bc_psi</name>
      <anchor>p65</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; RefineAlgorithm&lt; DIM &gt; &gt;</type>
      <name>d_fill_new_level</name>
      <anchor>p66</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; RefineAlgorithm&lt; DIM &gt; &gt;</type>
      <name>d_fill_bdry_compute_stable_dt</name>
      <anchor>p67</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Array&lt; Pointer&lt; RefineSchedule&lt; DIM &gt; &gt; &gt;</type>
      <name>d_fill_bdry_sched_compute_stable_dt</name>
      <anchor>p68</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Array&lt; Pointer&lt; RefineAlgorithm&lt; DIM &gt; &gt; &gt;</type>
      <name>d_fill_bdry_time_advance</name>
      <anchor>p69</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Array&lt; Array&lt; Pointer&lt; RefineSchedule&lt; DIM &gt; &gt; &gt; &gt;</type>
      <name>d_fill_bdry_sched_time_advance</name>
      <anchor>p70</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>LSMLIB::LevelSetFunctionIntegratorStrategy</name>
    <filename>classLSMLIB_1_1LevelSetFunctionIntegratorStrategy.html</filename>
    <templarg>DIM</templarg>
    <member kind="function" virtualness="pure">
      <type>virtual int</type>
      <name>getPhiPatchDataHandle</name>
      <anchor>z47_0</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual int</type>
      <name>getPsiPatchDataHandle</name>
      <anchor>z47_1</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual int</type>
      <name>getControlVolumePatchDataHandle</name>
      <anchor>z47_2</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual LSMLIB_REAL</type>
      <name>getStartTime</name>
      <anchor>z49_0</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual LSMLIB_REAL</type>
      <name>getEndTime</name>
      <anchor>z49_1</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual LSMLIB_REAL</type>
      <name>getCurrentTime</name>
      <anchor>z49_2</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>endTimeReached</name>
      <anchor>z49_3</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual int</type>
      <name>numIntegrationStepsTaken</name>
      <anchor>z49_4</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setBoundaryConditions</name>
      <anchor>z51_0</anchor>
      <arglist>(const  IntVector&lt; DIM &gt; &amp;lower_bc, const  IntVector&lt; DIM &gt; &amp;upper_bc, const  LEVEL_SET_FCN_TYPE level_set_fcn, const  int component=-1)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual LSMLIB_REAL</type>
      <name>computeStableDt</name>
      <anchor>z53_0</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>advanceLevelSetFunctions</name>
      <anchor>z53_1</anchor>
      <arglist>(const  LSMLIB_REAL dt)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual int</type>
      <name>getReinitializationInterval</name>
      <anchor>z55_0</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setReinitializationInterval</name>
      <anchor>z55_1</anchor>
      <arglist>(const  int interval)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>reinitializeLevelSetFunctions</name>
      <anchor>z55_2</anchor>
      <arglist>(const  LEVEL_SET_FCN_TYPE level_set_fcn=LSMLIB::PHI, const  int max_iterations=-1)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual int</type>
      <name>getOrthogonalizationInterval</name>
      <anchor>z55_3</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setOrthogonalizationInterval</name>
      <anchor>z55_4</anchor>
      <arglist>(const  int interval)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>orthogonalizeLevelSetFunctions</name>
      <anchor>z55_5</anchor>
      <arglist>(const  LEVEL_SET_FCN_TYPE level_set_fcn, const  int max_reinit_iterations=-1, const  int max_ortho_iterations=-1)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>initializeLevelData</name>
      <anchor>z57_0</anchor>
      <arglist>(const  Pointer&lt; BasePatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int level_number, const  double init_data_time, const  bool can_be_refined, const  bool initial_time, const  Pointer&lt; BasePatchLevel&lt; DIM &gt; &gt; old_level=Pointer&lt; BasePatchLevel&lt; DIM &gt; &gt;((0)), const  bool allocate_data=true)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>applyGradientDetector</name>
      <anchor>z57_1</anchor>
      <arglist>(const  Pointer&lt; BasePatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int level_number, const  double error_data_time, const  int tag_index, const  bool initial_time, const  bool uses_richardson_extrapolation_too)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>resetHierarchyConfiguration</name>
      <anchor>z57_2</anchor>
      <arglist>(Pointer&lt; BasePatchHierarchy&lt; DIM &gt; &gt; hierarchy, int coarsest_level, int finest_level)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>preprocessInitializeVelocityField</name>
      <anchor>z59_0</anchor>
      <arglist>(int &amp;phi_handle, int &amp;psi_handle, const  Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int level_number)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>postprocessInitializeVelocityField</name>
      <anchor>z59_1</anchor>
      <arglist>(const  Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int level_number)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LevelSetFunctionIntegratorStrategy</name>
      <anchor>z61_0</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>LSMLIB::LevelSetMethodAlgorithm</name>
    <filename>classLSMLIB_1_1LevelSetMethodAlgorithm.html</filename>
    <templarg>DIM</templarg>
    <member kind="function">
      <type></type>
      <name>LevelSetMethodAlgorithm</name>
      <anchor>z63_0</anchor>
      <arglist>(Pointer&lt; Database &gt; lsm_algorithm_input_db, Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; patch_hierarchy, LevelSetMethodPatchStrategy&lt; DIM &gt; *patch_strategy, LevelSetMethodVelocityFieldStrategy&lt; DIM &gt; *velocity_field_strategy, const  int num_level_set_fcn_components=1, const  int codimension=1, const  string &amp;object_name=&quot;LevelSetMethodAlgorithm&quot;)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LevelSetMethodAlgorithm</name>
      <anchor>z63_1</anchor>
      <arglist>(Pointer&lt; LevelSetFunctionIntegratorStrategy&lt; DIM &gt; &gt; lsm_integrator_strategy, Pointer&lt; LevelSetMethodGriddingStrategy&lt; DIM &gt; &gt; lsm_gridding_strategy, const  string &amp;object_name=&quot;LevelSetMethodAlgorithm&quot;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LevelSetMethodAlgorithm</name>
      <anchor>z63_2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getPhiPatchDataHandle</name>
      <anchor>z65_0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getPsiPatchDataHandle</name>
      <anchor>z65_1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getControlVolumePatchDataHandle</name>
      <anchor>z65_2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual LSMLIB_REAL</type>
      <name>getStartTime</name>
      <anchor>z67_0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual LSMLIB_REAL</type>
      <name>getEndTime</name>
      <anchor>z67_1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual LSMLIB_REAL</type>
      <name>getCurrentTime</name>
      <anchor>z67_2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>endTimeReached</name>
      <anchor>z67_3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>numIntegrationStepsTaken</name>
      <anchor>z67_4</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>printClassData</name>
      <anchor>z67_5</anchor>
      <arglist>(ostream &amp;os) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getSpatialDerivativeType</name>
      <anchor>z68_0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getSpatialDerivativeOrder</name>
      <anchor>z68_1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getTVDRungeKuttaOrder</name>
      <anchor>z68_2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual LSMLIB_REAL</type>
      <name>getCFLNumber</name>
      <anchor>z68_3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setBoundaryConditions</name>
      <anchor>z69_0</anchor>
      <arglist>(const  IntVector&lt; DIM &gt; &amp;lower_bc, const  IntVector&lt; DIM &gt; &amp;upper_bc, const  LEVEL_SET_FCN_TYPE level_set_fcn, const  int component=-1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>initializeLevelSetMethodCalculation</name>
      <anchor>z71_0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual LSMLIB_REAL</type>
      <name>computeStableDt</name>
      <anchor>z71_1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>advanceLevelSetFunctions</name>
      <anchor>z71_2</anchor>
      <arglist>(const  LSMLIB_REAL dt)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getReinitializationInterval</name>
      <anchor>z73_0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setReinitializationInterval</name>
      <anchor>z73_1</anchor>
      <arglist>(const  int interval)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>reinitializeLevelSetFunctions</name>
      <anchor>z73_2</anchor>
      <arglist>(const  LEVEL_SET_FCN_TYPE level_set_fcn=LSMLIB::PHI, const  int max_iterations=-1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getOrthogonalizationInterval</name>
      <anchor>z73_3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setOrthogonalizationInterval</name>
      <anchor>z73_4</anchor>
      <arglist>(const  int interval)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>orthogonalizeLevelSetFunctions</name>
      <anchor>z73_5</anchor>
      <arglist>(const  LEVEL_SET_FCN_TYPE level_set_fcn, const  int max_reinit_iterations=-1, const  int max_ortho_iterations=-1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetHierarchyConfiguration</name>
      <anchor>z75_0</anchor>
      <arglist>(const  Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int coarsest_level, const  int finest_level)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>regridPatchHierarchy</name>
      <anchor>z75_1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Pointer&lt; FieldExtensionAlgorithm&lt; DIM &gt; &gt;</type>
      <name>getFieldExtensionAlgorithm</name>
      <anchor>z77_0</anchor>
      <arglist>(Pointer&lt; Database &gt; input_db, const  int field_handle, const  LEVEL_SET_FCN_TYPE level_set_fcn, const  string &amp;object_name=&quot;FieldExtensionAlgorithm&quot;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Pointer&lt; FieldExtensionAlgorithm&lt; DIM &gt; &gt;</type>
      <name>getFieldExtensionAlgorithm</name>
      <anchor>z77_1</anchor>
      <arglist>(const  int field_handle, const  LEVEL_SET_FCN_TYPE level_set_fcn, SPATIAL_DERIVATIVE_TYPE spatial_derivative_type=UNKNOWN, int spatial_derivative_order=0, int tvd_runge_kutta_order=0, LSMLIB_REAL cfl_number=0, const  LSMLIB_REAL stop_distance=0.0, const  int max_iterations=0, const  LSMLIB_REAL iteration_stop_tolerance=0.0, const  bool verbose_mode=false, const  string &amp;object_name=&quot;FieldExtensionAlgorithm&quot;)</arglist>
    </member>
    <member kind="variable" protection="public">
      <type>string</type>
      <name>d_object_name</name>
      <anchor>p0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; LevelSetFunctionIntegratorStrategy&lt; DIM &gt; &gt;</type>
      <name>d_lsm_integrator_strategy</name>
      <anchor>p1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; LevelSetMethodGriddingStrategy&lt; DIM &gt; &gt;</type>
      <name>d_lsm_gridding_strategy</name>
      <anchor>p2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Array&lt; Pointer&lt; FieldExtensionAlgorithm&lt; DIM &gt; &gt; &gt;</type>
      <name>d_field_extension_alg_list</name>
      <anchor>p3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt;</type>
      <name>d_patch_hierarchy</name>
      <anchor>p4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_using_standard_level_set_fcn_integrator</name>
      <anchor>p5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>SPATIAL_DERIVATIVE_TYPE</type>
      <name>d_spatial_derivative_type</name>
      <anchor>p6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_spatial_derivative_order</name>
      <anchor>p7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_tvd_runge_kutta_order</name>
      <anchor>p8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LSMLIB_REAL</type>
      <name>d_cfl_number</name>
      <anchor>p9</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>LSMLIB::LevelSetMethodGriddingAlgorithm</name>
    <filename>classLSMLIB_1_1LevelSetMethodGriddingAlgorithm.html</filename>
    <templarg>DIM</templarg>
    <base>LSMLIB::LevelSetMethodGriddingStrategy</base>
    <member kind="function">
      <type></type>
      <name>LevelSetMethodGriddingAlgorithm</name>
      <anchor>z79_0</anchor>
      <arglist>(Pointer&lt; Database &gt; input_db, Pointer&lt; BasePatchHierarchy&lt; DIM &gt; &gt; patch_hierarchy, Pointer&lt; LevelSetFunctionIntegratorStrategy&lt; DIM &gt; &gt; lsm_integrator_strategy, const  string &amp;object_name=&quot;LevelSetMethodGriddingAlgorithm&quot;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LevelSetMethodGriddingAlgorithm</name>
      <anchor>z79_1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>registerVelocityFieldStrategy</name>
      <anchor>z81_0</anchor>
      <arglist>(LevelSetMethodVelocityFieldStrategy&lt; DIM &gt; *velocity_field_strategy)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>initializePatchHierarchy</name>
      <anchor>z83_0</anchor>
      <arglist>(const  LSMLIB_REAL time)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>regridPatchHierarchy</name>
      <anchor>z83_1</anchor>
      <arglist>(LSMLIB_REAL time)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>initializeLevelData</name>
      <anchor>z85_0</anchor>
      <arglist>(const  Pointer&lt; BasePatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int level_number, const  double init_data_time, const  bool can_be_refined, const  bool initial_time, const  Pointer&lt; BasePatchLevel&lt; DIM &gt; &gt; old_level=Pointer&lt; BasePatchLevel&lt; DIM &gt; &gt;(NULL), const  bool allocate_data=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>resetHierarchyConfiguration</name>
      <anchor>z85_1</anchor>
      <arglist>(const  Pointer&lt; BasePatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int coarsest_level, const  int finest_level)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>tagCellsForRefinement</name>
      <anchor>z85_2</anchor>
      <arglist>(const  Pointer&lt; BasePatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int level_number, const  double regrid_time, const  int tag_index, const  bool initial_time, const  bool coarsest_sync_level, const  bool can_be_refined, const  double regrid_start_time=0)</arglist>
    </member>
    <member kind="function" protection="public">
      <type>void</type>
      <name>getFromInput</name>
      <anchor>z87_0</anchor>
      <arglist>(Pointer&lt; Database &gt; input_db)</arglist>
    </member>
    <member kind="variable" protection="public">
      <type>string</type>
      <name>d_object_name</name>
      <anchor>p0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; BasePatchHierarchy&lt; DIM &gt; &gt;</type>
      <name>d_patch_hierarchy</name>
      <anchor>p1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; SAMRAI::mesh::GriddingAlgorithm&lt; DIM &gt; &gt;</type>
      <name>d_gridding_alg</name>
      <anchor>p2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_refine_boxes</name>
      <anchor>p3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_gradient_detector</name>
      <anchor>p4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_richardson_extrapolation</name>
      <anchor>p5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; LevelSetFunctionIntegratorStrategy&lt; DIM &gt; &gt;</type>
      <name>d_lsm_integrator_strategy</name>
      <anchor>p6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Array&lt; Pointer&lt; LevelSetMethodVelocityFieldStrategy&lt; DIM &gt; &gt; &gt;</type>
      <name>d_velocity_field_strategies</name>
      <anchor>p7</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>LSMLIB::LevelSetMethodGriddingStrategy</name>
    <filename>classLSMLIB_1_1LevelSetMethodGriddingStrategy.html</filename>
    <templarg>DIM</templarg>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>registerVelocityFieldStrategy</name>
      <anchor>z89_0</anchor>
      <arglist>(LevelSetMethodVelocityFieldStrategy&lt; DIM &gt; *velocity_field_strategy)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>initializePatchHierarchy</name>
      <anchor>z91_0</anchor>
      <arglist>(const  LSMLIB_REAL time)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>resetHierarchyConfiguration</name>
      <anchor>z91_1</anchor>
      <arglist>(const  Pointer&lt; BasePatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int coarsest_level, const  int finest_level)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>regridPatchHierarchy</name>
      <anchor>z91_2</anchor>
      <arglist>(LSMLIB_REAL time)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LevelSetMethodGriddingStrategy</name>
      <anchor>z93_0</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>LSMLIB::LevelSetMethodPatchStrategy</name>
    <filename>classLSMLIB_1_1LevelSetMethodPatchStrategy.html</filename>
    <templarg>DIM</templarg>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>initializeLevelSetFunctionsOnPatch</name>
      <anchor>z95_0</anchor>
      <arglist>(Patch&lt; DIM &gt; &amp;patch, const  LSMLIB_REAL time, const  int phi_handle, const  int psi_handle)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setLevelSetFunctionBoundaryConditions</name>
      <anchor>z95_1</anchor>
      <arglist>(Patch&lt; DIM &gt; &amp;patch, const  LSMLIB_REAL fill_time, const  int phi_handle, const  int psi_handle, const  IntVector&lt; DIM &gt; &amp;ghost_width_to_fill)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual LSMLIB_REAL</type>
      <name>computeStableDtOnPatch</name>
      <anchor>z97_0</anchor>
      <arglist>(Patch&lt; DIM &gt; &amp;patch, LevelSetFunctionIntegrator&lt; DIM &gt; *lsm_integrator, LevelSetMethodVelocityFieldStrategy&lt; DIM &gt; *velocity_field_strategy)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LevelSetMethodPatchStrategy</name>
      <anchor>z99_0</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>LSMLIB::LevelSetMethodToolbox</name>
    <filename>classLSMLIB_1_1LevelSetMethodToolbox.html</filename>
    <templarg>DIM</templarg>
    <member kind="enumeration">
      <name>UNIT_NORMAL_TYPE</name>
      <anchor>z101_0</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>PHI_UPWIND</name>
      <anchor>z101_0w0</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>AVERAGE</name>
      <anchor>z101_0w1</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" static="yes">
      <type>void</type>
      <name>computeUpwindSpatialDerivatives</name>
      <anchor>z103_0</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  SPATIAL_DERIVATIVE_TYPE spatial_derivative_type, const  int spatial_derivative_order, const  int grad_phi_handle, const  int phi_handle, const  int upwind_function_handle, const  int phi_component=0)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>void</type>
      <name>computePlusAndMinusSpatialDerivatives</name>
      <anchor>z103_1</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  SPATIAL_DERIVATIVE_TYPE spatial_derivative_type, const  int spatial_derivative_order, const  int grad_phi_plus_handle, const  int grad_phi_minus_handle, const  int phi_handle, const  int phi_component=0)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>void</type>
      <name>computeCentralSpatialDerivatives</name>
      <anchor>z103_2</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int spatial_derivative_order, const  int grad_phi_handle, const  int phi_handle, const  int phi_component=0)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>void</type>
      <name>TVDRK1Step</name>
      <anchor>z105_0</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int u_next_handle, const  int u_cur_handle, const  int rhs_handle, const  LSMLIB_REAL dt, const  int u_next_component=0, const  int u_cur_component=0, const  int rhs_component=0)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>void</type>
      <name>TVDRK2Stage1</name>
      <anchor>z105_1</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int u_stage1_handle, const  int u_cur_handle, const  int rhs_handle, const  LSMLIB_REAL dt, const  int u_stage1_component=0, const  int u_cur_component=0, const  int rhs_component=0)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>void</type>
      <name>TVDRK2Stage2</name>
      <anchor>z105_2</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int u_next_handle, const  int u_stage1_handle, const  int u_cur_handle, const  int rhs_handle, const  LSMLIB_REAL dt, const  int u_next_component=0, const  int u_stage1_component=0, const  int u_cur_component=0, const  int rhs_component=0)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>void</type>
      <name>TVDRK3Stage1</name>
      <anchor>z105_3</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int u_stage1_handle, const  int u_cur_handle, const  int rhs_handle, const  LSMLIB_REAL dt, const  int u_stage1_component=0, const  int u_cur_component=0, const  int rhs_component=0)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>void</type>
      <name>TVDRK3Stage2</name>
      <anchor>z105_4</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int u_stage2_handle, const  int u_stage1_handle, const  int u_cur_handle, const  int rhs_handle, const  LSMLIB_REAL dt, const  int u_stage2_component=0, const  int u_stage1_component=0, const  int u_cur_component=0, const  int rhs_component=0)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>void</type>
      <name>TVDRK3Stage3</name>
      <anchor>z105_5</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int u_next_handle, const  int u_stage2_handle, const  int u_cur_handle, const  int rhs_handle, const  LSMLIB_REAL dt, const  int u_next_component=0, const  int u_stage2_component=0, const  int u_cur_component=0, const  int rhs_component=0)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>void</type>
      <name>computeDistanceFunctionUsingFMM</name>
      <anchor>z107_0</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int spatial_derivative_order, const  int distance_function_handle, const  int phi_handle, const  int distance_function_component=0, const  int phi_component=0)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>void</type>
      <name>computeExtensionFieldsUsingFMM</name>
      <anchor>z107_1</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int spatial_derivative_order, const  vector&lt; int &gt; &amp;extension_field_handles, const  int distance_function_handle, const  vector&lt; int &gt; &amp;source_field_handles, const  int phi_handle, const  int extension_field_component=0, const  int distance_function_component=0, const  int source_field_component=0, const  int phi_component=0)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>void</type>
      <name>computeUnitNormalVectorFromPhi</name>
      <anchor>z109_0</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; patch_hierarchy, SPATIAL_DERIVATIVE_TYPE spatial_derivative_type, const  int spatial_derivative_order, const  int normal_vector_handle, const  int phi_handle, const  UNIT_NORMAL_TYPE unit_normal_type, const  int phi_component=0)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>void</type>
      <name>computeSignedUnitNormalVectorFromPhi</name>
      <anchor>z109_1</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; patch_hierarchy, SPATIAL_DERIVATIVE_TYPE spatial_derivative_type, const  int spatial_derivative_order, const  int normal_vector_handle, const  int phi_handle, const  UNIT_NORMAL_TYPE unit_normal_type, const  int phi_component=0)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>void</type>
      <name>computeUnitNormalVectorFromGradPhi</name>
      <anchor>z109_2</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; patch_hierarchy, const  int normal_vector_handle, const  int grad_phi_handle)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>void</type>
      <name>computeSignedUnitNormalVectorFromGradPhi</name>
      <anchor>z109_3</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; patch_hierarchy, const  int normal_vector_handle, const  int grad_phi_handle, const  int phi_handle, const  int phi_component=0)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>LSMLIB_REAL</type>
      <name>computeVolumeOfRegionDefinedByZeroLevelSet</name>
      <anchor>z109_4</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; patch_hierarchy, const  int phi_handle, const  int control_volume_handle, const  int region_indicator, const  int phi_component=0, const  int heaviside_width=3)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>LSMLIB_REAL</type>
      <name>computeVolumeOfZeroLevelSet</name>
      <anchor>z109_5</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; patch_hierarchy, const  int phi_handle, const  int grad_phi_handle, const  int control_volume_handle, const  int phi_component=0, const  int delta_width=3)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>LSMLIB_REAL</type>
      <name>computeVolumeIntegral</name>
      <anchor>z109_6</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; patch_hierarchy, const  int F_handle, const  int phi_handle, const  int control_volume_handle, const  int region_integrator, const  int F_component=0, const  int phi_component=0, const  int heaviside_width=3)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>LSMLIB_REAL</type>
      <name>computeSurfaceIntegral</name>
      <anchor>z109_7</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; patch_hierarchy, const  int F_handle, const  int phi_handle, const  int grad_phi_handle, const  int control_volume_handle, const  int F_component=0, const  int phi_component=0, const  int delta_width=3)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>LSMLIB_REAL</type>
      <name>computeStableAdvectionDt</name>
      <anchor>z111_0</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; patch_hierarchy, const  int velocity_handle, const  int control_volume_handle, const  LSMLIB_REAL cfl_number)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>LSMLIB_REAL</type>
      <name>computeStableNormalVelocityDt</name>
      <anchor>z111_1</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; patch_hierarchy, const  int normal_velocity_handle, const  int grad_phi_plus_handle, const  int grad_phi_minus_handle, const  int control_volume_handle, const  LSMLIB_REAL cfl_number)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>LSMLIB_REAL</type>
      <name>maxNormOfDifference</name>
      <anchor>z111_2</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int field1_handle, const  int field2_handle, const  int control_volume_handle, const  int field1_component=0, const  int field2_component=0)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>void</type>
      <name>computeControlVolumes</name>
      <anchor>z111_3</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; patch_hierarchy, int control_volume_handle)</arglist>
    </member>
    <member kind="function" static="yes">
      <type>void</type>
      <name>copySAMRAIData</name>
      <anchor>z111_4</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; patch_hierarchy, const  int dst_handle, const  int src_handle, const  int dst_component=0, const  int src_component=0)</arglist>
    </member>
    <member kind="function" protection="public" static="yes">
      <type>void</type>
      <name>initializeComputeSpatialDerivativesParameters</name>
      <anchor>z113_0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="public" static="yes">
      <type>void</type>
      <name>initializeComputeUnitNormalParameters</name>
      <anchor>z113_1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="public" static="yes">
      <type>int</type>
      <name>s_D1_one_ghostcell_handle</name>
      <anchor>t0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public" static="yes">
      <type>int</type>
      <name>s_D1_two_ghostcells_handle</name>
      <anchor>t1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public" static="yes">
      <type>int</type>
      <name>s_D2_two_ghostcells_handle</name>
      <anchor>t2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public" static="yes">
      <type>int</type>
      <name>s_D1_three_ghostcells_handle</name>
      <anchor>t3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public" static="yes">
      <type>int</type>
      <name>s_D2_three_ghostcells_handle</name>
      <anchor>t4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public" static="yes">
      <type>int</type>
      <name>s_D3_three_ghostcells_handle</name>
      <anchor>t5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public" static="yes">
      <type>int</type>
      <name>s_compute_normal_grad_phi_handle</name>
      <anchor>t6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public" static="yes">
      <type>int</type>
      <name>s_compute_normal_grad_phi_plus_handle</name>
      <anchor>t7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public" static="yes">
      <type>int</type>
      <name>s_compute_normal_grad_phi_minus_handle</name>
      <anchor>t8</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>LSMLIB::LevelSetMethodVelocityFieldStrategy</name>
    <filename>classLSMLIB_1_1LevelSetMethodVelocityFieldStrategy.html</filename>
    <templarg>DIM</templarg>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>providesExternalVelocityField</name>
      <anchor>z115_0</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>providesNormalVelocityField</name>
      <anchor>z115_1</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual int</type>
      <name>getExternalVelocityFieldPatchDataHandle</name>
      <anchor>z115_2</anchor>
      <arglist>(const  int component) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual int</type>
      <name>getNormalVelocityFieldPatchDataHandle</name>
      <anchor>z115_3</anchor>
      <arglist>(const  LEVEL_SET_FCN_TYPE level_set_fcn, const  int component) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>computeVelocityField</name>
      <anchor>z115_4</anchor>
      <arglist>(const  LSMLIB_REAL time, const  int phi_handle, const  int psi_handle, const  int component)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setCurrentTime</name>
      <anchor>z117_0</anchor>
      <arglist>(const  LSMLIB_REAL time)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual LSMLIB_REAL</type>
      <name>computeStableDt</name>
      <anchor>z117_1</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>initializeLevelData</name>
      <anchor>z119_0</anchor>
      <arglist>(const  Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int level_number, const  LSMLIB_REAL init_data_time, const  int phi_handle, const  int psi_handle, const  bool can_be_refined, const  bool initial_time, const  Pointer&lt; PatchLevel&lt; DIM &gt; &gt; old_level=Pointer&lt; PatchLevel&lt; DIM &gt; &gt;((0)), const  bool allocate_data=true)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>resetHierarchyConfiguration</name>
      <anchor>z121_0</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, int coarsest_level, int finest_level)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>tagCellsForRefinement</name>
      <anchor>z121_1</anchor>
      <arglist>(const  Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int level_number, const  int tag_handle)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LevelSetMethodVelocityFieldStrategy</name>
      <anchor>z123_0</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>LSMLIB::OrthogonalizationAlgorithm</name>
    <filename>classLSMLIB_1_1OrthogonalizationAlgorithm.html</filename>
    <templarg>DIM</templarg>
    <member kind="function">
      <type></type>
      <name>OrthogonalizationAlgorithm</name>
      <anchor>z125_0</anchor>
      <arglist>(Pointer&lt; Database &gt; input_db, Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int phi_handle, const  int psi_handle, const  int control_volume_handle, const  string &amp;object_name=&quot;OrthogonalizationAlgorithm&quot;, const  IntVector&lt; DIM &gt; &amp;phi_ghostcell_width=0, const  IntVector&lt; DIM &gt; &amp;psi_ghostcell_width=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>OrthogonalizationAlgorithm</name>
      <anchor>z125_1</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int phi_handle, const  int psi_handle, const  int control_volume_handle, const  SPATIAL_DERIVATIVE_TYPE spatial_derivative_type, const  int spatial_derivative_order, const  int tvd_runge_kutta_order, const  LSMLIB_REAL cfl_number, const  LSMLIB_REAL stop_distance=0.0, const  int max_iterations=0, const  LSMLIB_REAL iteration_stop_tolerance=0.0, const  bool verbose_mode=false, const  string &amp;object_name=&quot;OrthogonalizationAlgorithm&quot;, const  IntVector&lt; DIM &gt; &amp;phi_ghostcell_width=0, const  IntVector&lt; DIM &gt; &amp;psi_ghostcell_width=0)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~OrthogonalizationAlgorithm</name>
      <anchor>z125_2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>orthogonalizeLevelSetFunctions</name>
      <anchor>z127_0</anchor>
      <arglist>(const  LEVEL_SET_FCN_TYPE level_set_fcn, const  int max_iterations=-1, const  IntVector&lt; DIM &gt; &amp;lower_bc_fixed=IntVector&lt; DIM &gt;(-1), const  IntVector&lt; DIM &gt; &amp;upper_bc_fixed=IntVector&lt; DIM &gt;(-1), const  IntVector&lt; DIM &gt; &amp;lower_bc_evolved=IntVector&lt; DIM &gt;(-1), const  IntVector&lt; DIM &gt; &amp;upper_bc_evolved=IntVector&lt; DIM &gt;(-1))</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>orthogonalizeLevelSetFunctionForSingleComponent</name>
      <anchor>z127_1</anchor>
      <arglist>(const  LEVEL_SET_FCN_TYPE level_set_fcn, const  int component=0, const  int max_iterations=-1, const  IntVector&lt; DIM &gt; &amp;lower_bc_fixed=IntVector&lt; DIM &gt;(-1), const  IntVector&lt; DIM &gt; &amp;upper_bc_fixed=IntVector&lt; DIM &gt;(-1), const  IntVector&lt; DIM &gt; &amp;lower_bc_evolved=IntVector&lt; DIM &gt;(-1), const  IntVector&lt; DIM &gt; &amp;upper_bc_evolved=IntVector&lt; DIM &gt;(-1))</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>resetHierarchyConfiguration</name>
      <anchor>z129_0</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int coarsest_level, const  int finest_level)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>getFromInput</name>
      <anchor>z131_0</anchor>
      <arglist>(Pointer&lt; Database &gt; db)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>checkParameters</name>
      <anchor>z131_1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="public">
      <type>string</type>
      <name>d_object_name</name>
      <anchor>p0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>SPATIAL_DERIVATIVE_TYPE</type>
      <name>d_spatial_derivative_type</name>
      <anchor>p1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_spatial_derivative_order</name>
      <anchor>p2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_tvd_runge_kutta_order</name>
      <anchor>p3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LSMLIB_REAL</type>
      <name>d_cfl_number</name>
      <anchor>p4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LSMLIB_REAL</type>
      <name>d_stop_distance</name>
      <anchor>p5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_max_iterations</name>
      <anchor>p6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LSMLIB_REAL</type>
      <name>d_iteration_stop_tol</name>
      <anchor>p7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_verbose_mode</name>
      <anchor>p8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt;</type>
      <name>d_patch_hierarchy</name>
      <anchor>p9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; CartesianGridGeometry&lt; DIM &gt; &gt;</type>
      <name>d_grid_geometry</name>
      <anchor>p10</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; FieldExtensionAlgorithm&lt; DIM &gt; &gt;</type>
      <name>d_fixed_phi_field_ext_alg</name>
      <anchor>p11</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; FieldExtensionAlgorithm&lt; DIM &gt; &gt;</type>
      <name>d_fixed_psi_field_ext_alg</name>
      <anchor>p12</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_phi_handle</name>
      <anchor>p13</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_psi_handle</name>
      <anchor>p14</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_control_volume_handle</name>
      <anchor>p15</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_stop_distance</name>
      <anchor>p16</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_max_iterations</name>
      <anchor>p17</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_iteration_stop_tol</name>
      <anchor>p18</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_num_field_components</name>
      <anchor>p19</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>LSMLIB::ReinitializationAlgorithm</name>
    <filename>classLSMLIB_1_1ReinitializationAlgorithm.html</filename>
    <templarg>DIM</templarg>
    <member kind="function">
      <type></type>
      <name>ReinitializationAlgorithm</name>
      <anchor>z133_0</anchor>
      <arglist>(Pointer&lt; Database &gt; input_db, Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int phi_handle, const  int control_volume_handle, const  string &amp;object_name=&quot;ReinitializationAlgorithm&quot;)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ReinitializationAlgorithm</name>
      <anchor>z133_1</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int phi_handle, const  int control_volume_handle, const  SPATIAL_DERIVATIVE_TYPE spatial_derivative_type, const  int spatial_derivative_order, const  int tvd_runge_kutta_order, const  LSMLIB_REAL cfl_number, const  LSMLIB_REAL stop_distance=0.0, const  int max_iterations=0, const  LSMLIB_REAL iteration_stop_tolerance=0.0, const  bool verbose_mode=false, const  string &amp;object_name=&quot;ReinitializationAlgorithm&quot;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~ReinitializationAlgorithm</name>
      <anchor>z133_2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>reinitializeLevelSetFunctions</name>
      <anchor>z135_0</anchor>
      <arglist>(const  int max_iterations=-1, const  IntVector&lt; DIM &gt; &amp;lower_bc=IntVector&lt; DIM &gt;(-1), const  IntVector&lt; DIM &gt; &amp;upper_bc=IntVector&lt; DIM &gt;(-1))</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>reinitializeLevelSetFunctionForSingleComponent</name>
      <anchor>z135_1</anchor>
      <arglist>(const  int component=0, const  int max_iterations=-1, const  IntVector&lt; DIM &gt; &amp;lower_bc=IntVector&lt; DIM &gt;(-1), const  IntVector&lt; DIM &gt; &amp;upper_bc=IntVector&lt; DIM &gt;(-1))</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>resetHierarchyConfiguration</name>
      <anchor>z137_0</anchor>
      <arglist>(Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt; hierarchy, const  int coarsest_level, const  int finest_level)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>advanceReinitializationEqnUsingTVDRK1</name>
      <anchor>z139_0</anchor>
      <arglist>(const  LSMLIB_REAL dt, const  int component, const  IntVector&lt; DIM &gt; &amp;lower_bc, const  IntVector&lt; DIM &gt; &amp;upper_bc)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>advanceReinitializationEqnUsingTVDRK2</name>
      <anchor>z139_1</anchor>
      <arglist>(const  LSMLIB_REAL dt, const  int component, const  IntVector&lt; DIM &gt; &amp;lower_bc, const  IntVector&lt; DIM &gt; &amp;upper_bc)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>advanceReinitializationEqnUsingTVDRK3</name>
      <anchor>z139_2</anchor>
      <arglist>(const  LSMLIB_REAL dt, const  int component, const  IntVector&lt; DIM &gt; &amp;lower_bc, const  IntVector&lt; DIM &gt; &amp;upper_bc)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>computeReinitializationEqnRHS</name>
      <anchor>z139_3</anchor>
      <arglist>(const  int phi_handle)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>initializeVariables</name>
      <anchor>z141_0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>initializeCommunicationObjects</name>
      <anchor>z141_1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>getFromInput</name>
      <anchor>z141_2</anchor>
      <arglist>(Pointer&lt; Database &gt; db)</arglist>
    </member>
    <member kind="function" protection="public" virtualness="virtual">
      <type>virtual void</type>
      <name>checkParameters</name>
      <anchor>z141_3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="public">
      <type>string</type>
      <name>d_object_name</name>
      <anchor>p0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>SPATIAL_DERIVATIVE_TYPE</type>
      <name>d_spatial_derivative_type</name>
      <anchor>p1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_spatial_derivative_order</name>
      <anchor>p2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_tvd_runge_kutta_order</name>
      <anchor>p3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LSMLIB_REAL</type>
      <name>d_cfl_number</name>
      <anchor>p4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LSMLIB_REAL</type>
      <name>d_stop_distance</name>
      <anchor>p5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_max_iterations</name>
      <anchor>p6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>LSMLIB_REAL</type>
      <name>d_iteration_stop_tol</name>
      <anchor>p7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_verbose_mode</name>
      <anchor>p8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; PatchHierarchy&lt; DIM &gt; &gt;</type>
      <name>d_patch_hierarchy</name>
      <anchor>p9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; CartesianGridGeometry&lt; DIM &gt; &gt;</type>
      <name>d_grid_geometry</name>
      <anchor>p10</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_phi_handle</name>
      <anchor>p11</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_control_volume_handle</name>
      <anchor>p12</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>vector&lt; int &gt;</type>
      <name>d_phi_scr_handles</name>
      <anchor>p13</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>IntVector&lt; DIM &gt;</type>
      <name>d_phi_scratch_ghostcell_width</name>
      <anchor>p14</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_rhs_handle</name>
      <anchor>p15</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_grad_phi_plus_handle</name>
      <anchor>p16</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_grad_phi_minus_handle</name>
      <anchor>p17</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_stop_distance</name>
      <anchor>p18</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_max_iterations</name>
      <anchor>p19</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_use_iteration_stop_tol</name>
      <anchor>p20</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>bool</type>
      <name>d_hierarchy_configuration_needs_reset</name>
      <anchor>p21</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>int</type>
      <name>d_num_phi_components</name>
      <anchor>p22</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>ComponentSelector</type>
      <name>d_scratch_data</name>
      <anchor>p23</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Pointer&lt; BoundaryConditionModule&lt; DIM &gt; &gt;</type>
      <name>d_bc_module</name>
      <anchor>p24</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Array&lt; Pointer&lt; RefineAlgorithm&lt; DIM &gt; &gt; &gt;</type>
      <name>d_phi_fill_bdry_alg</name>
      <anchor>p25</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="public">
      <type>Array&lt; Array&lt; Pointer&lt; RefineSchedule&lt; DIM &gt; &gt; &gt; &gt;</type>
      <name>d_phi_fill_bdry_sched</name>
      <anchor>p26</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>mesh</name>
    <filename>namespacemesh.html</filename>
  </compound>
  <compound kind="namespace">
    <name>SAMRAI</name>
    <filename>namespaceSAMRAI.html</filename>
  </compound>
  <compound kind="namespace">
    <name>std</name>
    <filename>namespacestd.html</filename>
  </compound>
  <compound kind="namespace">
    <name>tbox</name>
    <filename>namespacetbox.html</filename>
  </compound>
  <compound kind="namespace">
    <name>xfer</name>
    <filename>namespacexfer.html</filename>
  </compound>
</tagfile>
