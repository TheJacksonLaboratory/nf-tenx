

class WorkflowMain {

    public static String logo(workflow) {
        """\n
  
        ------------------------------------------------------
  
        ███╗   ██╗███████╗ ████████╗███████╗███╗   ██╗██╗  ██╗
        ████╗  ██║██╔════╝ ╚══██╔══╝██╔════╝████╗  ██║╚██╗██╔╝
        ██╔██╗ ██║█████╗█████╗██║   █████╗  ██╔██╗ ██║ ╚███╔╝ 
        ██║╚██╗██║██╔══╝╚════╝██║   ██╔══╝  ██║╚██╗██║ ██╔██╗ 
        ██║ ╚████║██║         ██║   ███████╗██║ ╚████║██╔╝ ██╗
        ╚═╝  ╚═══╝╚═╝         ╚═╝   ╚══════╝╚═╝  ╚═══╝╚═╝  ╚═╝
  
        $workflow.manifest.name  v$workflow.manifest.version
        ------------------------------------------------------
  
        """.stripIndent()
  
    }
  
    public static String help(workflow, params, log) {
        def help_string = ""
        help_string += logo(workflow)
        return help_string
    }
  
    public static void initialize(workflow, params, log) {
        if (params.help) {
            log.info help(workflow, params, log)
            System.exit(0)
        }

        log.info logo(workflow)

        if (!params.samplesheet) {
            log.error "Please specify a samplesheet for the pipeline with '--samplesheet path/to/samplesheet.yml'"
            System.exit(1)
        }
        
    }
}
