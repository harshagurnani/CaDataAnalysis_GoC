macro crop_stack2{ 
makeRectangle(1,1,2,2);
run("Specify...", "width=324 height=211 x=260 y=34 slice=1");
run("Duplicate...", "title=whiskpadR_wcam duplicate")
}