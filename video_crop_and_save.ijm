macro split_stack{ 



makeRectangle(0,680,358,218);
run("Specify...", "width=109 height=125 x=395 y=509 slice=1");
run("Duplicate...", "title=whiskpadR_ecam duplicate")

makeRectangle(1,1,2,2);
run("Specify...", "width=358 height=218 x=0 y=680 slice=1");
run("Duplicate...", "title=forepawR_ecam duplicate")
}