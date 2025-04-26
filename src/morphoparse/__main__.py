from morphoparse import cli, core

def main():
    args = cli.get_args()
    core.run_pipeline(args)

if __name__ == "__main__":
    main()