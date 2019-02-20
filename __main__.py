# this file is only called when the package is called from the command line
import argparse
from .alginment import smith_waterman

def get_args():
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("-t","--total",type=int,default="100",help='total number of lists to have')
	return parser.parse_args()





return start,end

def main():

	args = get_args()
    starting_position, length_string = smith_waterman(a, b)

if __name__ == "__main__":
	main()
