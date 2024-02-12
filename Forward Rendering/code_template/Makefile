rasterizer:
	g++ *.cpp -std=c++11 -O3 -o rasterizer

clipping: rasterizer
	./rasterizer ../clipping_example/empty_box_clipped.xml
	mv *.ppm ../clipping_example/rasterized_outputs/

culling_disabled: rasterizer
	./rasterizer ../culling_disabled_inputs/empty_box.xml
	./rasterizer ../culling_disabled_inputs/filled_box.xml
	./rasterizer ../culling_disabled_inputs/flag_brazil.xml
	./rasterizer ../culling_disabled_inputs/flag_czechia.xml
	./rasterizer ../culling_disabled_inputs/flag_czechia_alternative.xml
	./rasterizer ../culling_disabled_inputs/flag_eu.xml
	./rasterizer ../culling_disabled_inputs/flag_eu_alternative.xml
	./rasterizer ../culling_disabled_inputs/flag_germany.xml
	./rasterizer ../culling_disabled_inputs/flag_turkey.xml
	./rasterizer ../culling_disabled_inputs/flag_turkey_alternative.xml
	./rasterizer ../culling_disabled_inputs/horse_and_mug.xml
	./rasterizer ../culling_disabled_inputs/sample.xml
	mv *.ppm ../culling_disabled_inputs/rasterized_outputs/

culling_enabled: rasterizer
	./rasterizer ../culling_enabled_inputs/empty_box.xml
	./rasterizer ../culling_enabled_inputs/filled_box.xml
	./rasterizer ../culling_enabled_inputs/flag_brazil.xml
	./rasterizer ../culling_enabled_inputs/flag_czechia.xml
	./rasterizer ../culling_enabled_inputs/flag_czechia_alternative.xml
	./rasterizer ../culling_enabled_inputs/flag_eu.xml
	./rasterizer ../culling_enabled_inputs/flag_eu_alternative.xml
	./rasterizer ../culling_enabled_inputs/flag_germany.xml
	./rasterizer ../culling_enabled_inputs/flag_iceland.xml
	./rasterizer ../culling_enabled_inputs/flag_turkey.xml
	./rasterizer ../culling_enabled_inputs/flag_turkey_alternative.xml
	./rasterizer ../culling_enabled_inputs/horse_and_mug.xml
	./rasterizer ../culling_enabled_inputs/sample.xml
	mv *.ppm ../culling_enabled_inputs/rasterized_outputs/

different_projection: rasterizer
	./rasterizer ../different_projection_type/flag_turkey/flag_turkey_orthographic.xml
	./rasterizer ../different_projection_type/flag_turkey/flag_turkey_perspective.xml
	./rasterizer ../different_projection_type/horse_and_mug/horse_and_mug_orthographic.xml
	./rasterizer ../different_projection_type/horse_and_mug/horse_and_mug_perspective.xml
	mv *.ppm ../different_projection_type/rasterized_outputs/
clean:
	rm -f rasterizer

clear:
	rm *.ppm