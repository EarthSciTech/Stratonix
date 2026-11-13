Stratonix export
================

Model:              synthetic_fold_two_faults v1.0.0 - Synthetic asymmetric anticline with a blind and an upward-reaching fault
Evidence basis:     DESIGN (design model, not an interpretation)
                    Every parameter of the description carries a value, unit, evidence class and source (see manifest.json).
Coordinate system:  local engineering frame
Vertical datum:     model datum, z = 0
Units:              metres (x, y, z, MD)
Z convention:       elevation, negative below the datum, in every file (TSurf: ZPOSITIVE Elevation)
Grids:              origin (0.00, 0.00), spacing 250 x 250 m, 33 columns x 25 rows, no rotation
Undefined values:   IRAP 9999900.0, ZMAP+ -99999.0 (fault gaps and outside the horizon)
Line separator:     IRAP lines use '999.000000 999.000000 999.000000' between lines
Fault sticks:       every 250 m along strike
Written by:         Stratonix 0.8.0

Files
-----
horizons/H1_faulted.irap                       H1 faulted surface; fault gaps undefined (9999900.0)
                                               format: IRAP classic ASCII grid; Petrel import: IRAP classic grid (ASCII)
horizons/H2_faulted.irap                       H2 faulted surface; fault gaps undefined (9999900.0)
                                               format: IRAP classic ASCII grid; Petrel import: IRAP classic grid (ASCII)
horizons/per_block/H2_block1.irap              H2 in block 1 (F1-HW_F2-FW), extended 2 nodes across the gap
                                               format: IRAP classic ASCII grid; Petrel import: IRAP classic grid (ASCII)
horizons/per_block/H2_block2.irap              H2 in block 2 (F1-FW_F2-HW), extended 2 nodes across the gap
                                               format: IRAP classic ASCII grid; Petrel import: IRAP classic grid (ASCII)
horizons/per_block/H2_block3.irap              H2 in block 3 (F1-HW_F2-HW), extended 2 nodes across the gap
                                               format: IRAP classic ASCII grid; Petrel import: IRAP classic grid (ASCII)
horizons/H3_faulted.irap                       H3 faulted surface; fault gaps undefined (9999900.0)
                                               format: IRAP classic ASCII grid; Petrel import: IRAP classic grid (ASCII)
horizons/per_block/H3_block1.irap              H3 in block 1 (F1-HW_F2-FW), extended 2 nodes across the gap
                                               format: IRAP classic ASCII grid; Petrel import: IRAP classic grid (ASCII)
horizons/per_block/H3_block2.irap              H3 in block 2 (F1-FW_F2-HW), extended 2 nodes across the gap
                                               format: IRAP classic ASCII grid; Petrel import: IRAP classic grid (ASCII)
horizons/per_block/H3_block3.irap              H3 in block 3 (F1-HW_F2-HW), extended 2 nodes across the gap
                                               format: IRAP classic ASCII grid; Petrel import: IRAP classic grid (ASCII)
horizons/H4_faulted.irap                       H4 faulted surface; fault gaps undefined (9999900.0)
                                               format: IRAP classic ASCII grid; Petrel import: IRAP classic grid (ASCII)
horizons/per_block/H4_block1.irap              H4 in block 1 (F1-HW_F2-FW), extended 2 nodes across the gap
                                               format: IRAP classic ASCII grid; Petrel import: IRAP classic grid (ASCII)
horizons/per_block/H4_block2.irap              H4 in block 2 (F1-FW_F2-HW), extended 2 nodes across the gap
                                               format: IRAP classic ASCII grid; Petrel import: IRAP classic grid (ASCII)
horizons/per_block/H4_block3.irap              H4 in block 3 (F1-HW_F2-HW), extended 2 nodes across the gap
                                               format: IRAP classic ASCII grid; Petrel import: IRAP classic grid (ASCII)
horizons/H5_faulted.irap                       H5 faulted surface; fault gaps undefined (9999900.0)
                                               format: IRAP classic ASCII grid; Petrel import: IRAP classic grid (ASCII)
horizons/H1_faulted.zmap                       H1 faulted surface; fault gaps null (-99999.0)
                                               format: ZMAP+ grid; Petrel import: Zmap+ grid (ASCII)
horizons/H2_faulted.zmap                       H2 faulted surface; fault gaps null (-99999.0)
                                               format: ZMAP+ grid; Petrel import: Zmap+ grid (ASCII)
horizons/H3_faulted.zmap                       H3 faulted surface; fault gaps null (-99999.0)
                                               format: ZMAP+ grid; Petrel import: Zmap+ grid (ASCII)
horizons/H4_faulted.zmap                       H4 faulted surface; fault gaps null (-99999.0)
                                               format: ZMAP+ grid; Petrel import: Zmap+ grid (ASCII)
horizons/H5_faulted.zmap                       H5 faulted surface; fault gaps null (-99999.0)
                                               format: ZMAP+ grid; Petrel import: Zmap+ grid (ASCII)
horizons/fault_polygons/H2_fault_polygons.dat  H2 fault polygons (F2): footwall cut-off then hanging-wall cut-off, closed
                                               format: IRAP classic lines; Petrel import: Irap classic lines (ASCII), as fault polygons
horizons/fault_polygons/H3_fault_polygons.dat  H3 fault polygons (F2): footwall cut-off then hanging-wall cut-off, closed
                                               format: IRAP classic lines; Petrel import: Irap classic lines (ASCII), as fault polygons
horizons/fault_polygons/H4_fault_polygons.dat  H4 fault polygons (F1, F2): footwall cut-off then hanging-wall cut-off, closed
                                               format: IRAP classic lines; Petrel import: Irap classic lines (ASCII), as fault polygons
faults/F1_sticks.dat                           F1: 14 sticks at 250 m
                                               format: IRAP classic lines; Petrel import: Irap classic lines (ASCII), as fault sticks
faults/F2_sticks.dat                           F2: 24 sticks at 250 m
                                               format: IRAP classic lines; Petrel import: Irap classic lines (ASCII), as fault sticks
faults/fault_sticks_points.txt                 fault stick vertices of all faults, attributes Fault and Stick
                                               format: Petrel points with attributes; Petrel import: Petrel points with attributes (ASCII)
faults/F1.ts                                   F1 triangulated plane clipped to the tip line (221 triangles)
                                               format: GOCAD TSurf; Petrel import: GOCAD TSurf (ASCII)
faults/F2.ts                                   F2 triangulated plane clipped to the tip line (559 triangles)
                                               format: GOCAD TSurf; Petrel import: GOCAD TSurf (ASCII)
horizons/H1_faulted.ts                         H1 faulted surface as triangles (no triangles across gaps)
                                               format: GOCAD TSurf; Petrel import: GOCAD TSurf (ASCII)
horizons/H2_faulted.ts                         H2 faulted surface as triangles (no triangles across gaps)
                                               format: GOCAD TSurf; Petrel import: GOCAD TSurf (ASCII)
horizons/H3_faulted.ts                         H3 faulted surface as triangles (no triangles across gaps)
                                               format: GOCAD TSurf; Petrel import: GOCAD TSurf (ASCII)
horizons/H4_faulted.ts                         H4 faulted surface as triangles (no triangles across gaps)
                                               format: GOCAD TSurf; Petrel import: GOCAD TSurf (ASCII)
horizons/H5_faulted.ts                         H5 faulted surface as triangles (no triangles across gaps)
                                               format: GOCAD TSurf; Petrel import: GOCAD TSurf (ASCII)
wells/well_heads.txt                           3 wells, KB elevation, TD as MD
                                               format: Petrel well heads; Petrel import: Well heads
wells/well_tops.txt                            15 model tops (horizon crossings) on all wells
                                               format: Petrel well tops; Petrel import: Petrel well tops (ASCII)

Fault blocks (per-block surfaces are named <horizon>_block<code>)
----------------------------------------------------------------
block1: F1-HW_F2-FW
block2: F1-FW_F2-HW
block3: F1-HW_F2-HW

The Petrel import options above are the file types to choose in Petrel's import dialog. Their behaviour
in Petrel is to be confirmed with the Petrel import checklist (STX-TN-02).
