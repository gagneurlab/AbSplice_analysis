workdir:
    "../../../"
    
configfile: "general_workflow/config_als.yaml"
    
import yaml
config_static_path = "general_workflow/config_static_structure.yaml"
with open(config_static_path, "r") as fd:
    config_static = yaml.safe_load(fd)
config_precomputed_path = "general_workflow/config_provided_files.yaml"
with open(config_precomputed_path, "r") as fd:
    config_precomputed = yaml.safe_load(fd)
    
OUTPUT_DIR_JUNCTION_ANNO = config['output_dir']['junction_annotation']
OUTPUT_DIR_BENCHMARK = config['output_dir']['benchmark']
OUTPUT_DIR_SPLICING = config['output_dir']['splicing_preds']
OUTPUT_DIR_OUTLIER = config['output_dir']['outtliers']
OUTPUT_DIR_VCF_ANNO = config['output_dir']['vcf_annotation']


rule als_compare_splicemap_motor_neuron_with_gtex_brain:
    input:
        als_splicemap_psi5 = config['splicemap']['psi5'].format(
            method=config['method'], event_filter=config['event_filter']
        ),
        als_splicemap_psi3 = config['splicemap']['psi3'].format(
            method=config['method'], event_filter=config['event_filter']
        ),
        gtex_splicemap_5 = expand(config_precomputed['splicemap_gtex']['psi5'],
            genome=config['genome'], gtex_version=config['gtex_version'], tissue=config['gtex_tissues'], 
            method=config['method'], event_filter=config['event_filter']
        ),
        gtex_splicemap_3 = expand(config_precomputed['splicemap_gtex']['psi3'],
            genome=config['genome'], gtex_version=config['gtex_version'], tissue=config['gtex_tissues'], 
            method=config['method'], event_filter=config['event_filter']
        ),
    output:
        shared_splice_sites = OUTPUT_DIR_BENCHMARK + config['gtex_comparison']['compare_splicemap']['shared_splice_sites'],
        ref_psi_correlation_psi5 = OUTPUT_DIR_BENCHMARK + config['gtex_comparison']['compare_splicemap']['ref_psi_correlation']['psi5'],
        ref_psi_correlation_psi3 = OUTPUT_DIR_BENCHMARK + config['gtex_comparison']['compare_splicemap']['ref_psi_correlation']['psi3'],
    notebook:
        "./compare_splicemaps.ipynb"     

        
