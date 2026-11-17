#coding=utf-8
from __future__ import print_function
import argparse
example_text = '''example:
    python this.py 
'''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]', epilog=example_text)

parser.add_argument('--sample_id','-s',type=str,help="sample id",required= True,metavar='')
parser.add_argument('--output', '-o',type= str,help="the output file",required= True,metavar='')
parser.add_argument('--input', '-i',type= str,help="the input template of config file",required= True,metavar='')
parser.add_argument('--path', '-p',type= str,help="the paht of the output file",required= True,metavar='')

def is_last_character_letter(input_string):
    """
    Check if the last character of the input string is a letter or not.
    and return the input string and the corresponding read length and depth.
    """
    if input_string[-1].isdigit():
        return input_string, 100, 7
    elif input_string[-1].isalpha():
        return input_string, 150, 30
    else:
        print("the sampleID is ERROR, The last character of the input string should be a digit or a letter.")

if __name__ == '__main__':
    args = parser.parse_args()
    sample_id=args.sample_id
    output=args.output
    input=args.input
    path=args.path
    id, read_length, depth = is_last_character_letter(sample_id)

    with open(output, 'w') as outf:
        with open(input, 'r') as configf:
            for line in configf:
                if line.strip(" ").startswith("id"):
                    print('    id: "'+sample_id+'"',file = outf)
                elif line.strip(" ").startswith("depth"):
                    print('    depth: '+str(depth),file = outf)
                elif line.strip(" ").startswith("read_length"):
                    print('    read_length: '+str(read_length),file = outf)
                elif line.strip(" ").startswith("PWD"):
                    print('    PWD: '+str(path),file = outf)
                elif line.strip(" ").startswith("pwd"):
                    print('pwd='+str(path),file = outf)
                else:
                    outf.write(line)