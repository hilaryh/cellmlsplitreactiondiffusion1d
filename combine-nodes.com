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
	

# gfx edit scene
gfx create window 1
# gfx edit spectrum
# gfx list node 5
gfx modify g_element "/" general clear;
gfx modify g_element "/" points domain_nodes coordinate coordinates tessellation default_points LOCAL glyph sphere size "0.1*0.1*0.1" offset 0,0,0 font default select_on material default data dependent spectrum default selected_material default_selected render_shaded;
gfx mod win 1 layout width 500 height 500
gfx modify window 1 image view_all
gfx modify spectrum default linear reverse range 0 0.00055 extend_above extend_below rainbow colour_range 0 1 component 1;


# gfx timekeeper default play speed 1 skip;
# gfx create time_editor

for($t=0;$t<$numfiles;$t++){
    gfx timekeeper default set_time $t
        $image = "TIME_STEP_SPEC_1.".sprintf("%04d",$t).".jpg"
        cmiss("gfx print jpg file $image")
    }
