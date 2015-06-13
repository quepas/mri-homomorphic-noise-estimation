from configuration import build_default
from experiment import run_example_experiments
import argparse

__author__ = 'quepas, jakubsieradzki'

def run_program():
    parser = argparse.ArgumentParser(description = 'Run MRI homomorphic noise estimation.',
                                     epilog = 'Coders: quepas, jakubsieradzki')
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-e", "--example",
                       help = "Run four predefined examples.",
                       action = 'store_true')
    group.add_argument("-c", "--config",
                       help = "Run experiment using configuration file.",
                       action = 'store')
    args = parser.parse_args()
    if args.example:
        config = build_default()
        run_example_experiments(config)
    elif args.config:
        print "Runing config"
    else:
        print "Use help option [-h] to get some help."

if __name__ == '__main__':
    run_program()
