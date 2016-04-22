opendir(datastuff,"/Users/hhunt1/Documents/sims/opencmiss/classicalfield_reactiondiffusion_cellmlsplitreactiondiffusion1d/output")
$nodefolder = "/Users/hhunt1/Documents/sims/opencmiss/classicalfield_reactiondiffusion_cellmlsplitreactiondiffusion1d/output/";

while ($name = readdir(datastuff)) {
	if (index($name,".exnode")!=-1 and index($name,"TIME_STEP_SPEC_1") !=-1) {
	    print $name
		print "\n"
		push(@exnodes,$name);
		
		}
	}
close(datastuff);
@exnodes[0] = $nodefolder . @exnodes[0];
cmiss("gfx read node @exnodes[0] time 0");
$numfiles=scalar @exnodes;
for ($i=1;$i<$numfiles;$i++) {
	@exnodes[$i] = $nodefolder . @exnodes[$i];
	cmiss("gfx read node @exnodes[$i] time $i");
	}
	
print "@exnodes"
gfx create window
gfx edit scene
gfx edit spectrum
gfx list node 5
gfx modify g_element "/" general clear;
gfx modify g_element "/" points domain_nodes coordinate coordinates tessellation default_points LOCAL glyph sphere size "1*1*1" offset 0,0,0 font default select_on material default data dependent spectrum default selected_material default_selected render_shaded;
