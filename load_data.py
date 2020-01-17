#!/usr/bin/env python
import pandas as pd
import requests
import bs4

# Basic data loader functions for a dataframe from a file
#    using some or all columns.
# Designed specifically for MaxQuant's proteinGroups
#More information may be found at:
#    https://github.com/PayneLab/SingleCellTMTQualityControl


def download_file(download_to_path="data/datafile", url_file_path="data/url.txt", password_file_path="data/password.txt"):
    """Download a file from a given url to the specified location.

    Parameters:
    path (str): The path to the file to save the file to on the local machine.

    Returns:
    str: The path the file was downloaded to.
    """

    #import pdb; pdb.set_trace()
    url_file = open(url_file_path, 'r')
    url = url_file.read().strip()
    url_file.close()

    password_file = open(password_file_path, 'r')
    password = password_file.read().strip()
    password_file.close()



    for i in range(2):

        with requests.Session() as session: # Use a session object to save cookies
            # Construct the urls for our GET and POST requests
            get_url = url
            post_url = get_url.replace("https://byu.box.com/shared", "https://byu.app.box.com/public")

            # Send initial GET request and parse the request token out of the response
            get_response = session.get(get_url)
            soup = bs4.BeautifulSoup(get_response.text, "html.parser")
            token_tag = soup.find(id="request_token")
            token = token_tag.get("value")

            # Send a POST request, with the password and token, to get the data
            payload = {
                'password': password,
                'request_token': token}
            response = session.post(post_url, data=payload)

            with open(download_to_path, 'wb') as dest:
                dest.write(response.content)
    return download_to_path


from tkinter import filedialog, Tk
def find_file():
    #Opens the dialog to select a file
    #returns the file path.
    r=Tk()
    r.withdraw()
    input_file =  filedialog.askopenfilename(initialdir = "/",title = "Select file", \
                                             filetypes = (("txt files","*.txt"),("all files","*.*")))
    r.destroy()
    return input_file

def print_columns(file):
    #Displays all column names for the user.
    #This is not used in plotting,
    # but provided as a convienience.
    with open(file, 'r') as _file:
        line = _file.readline().strip()
        headings = line.split('\t')
    print (headings)
    return headings

def get_cols(file, prefix=None, experiment=None, contains=[], not_contains=[]):
    #Generates a list of column names fitting certain requirements.
    #    prefix is what they must start with.
    #    experiment is what they end with
    #    contains is a list of phrases that
    #        must be somewhere in the name
    #    not_contains is a list of phrases
    #        that must not be in the name.
    #This function is used by load_dataset()
    with open(file, 'r') as _file:
        line = _file.readline().strip()
        headings = line.split('\t')
    headings = [i.strip('"') for i in headings]
    if prefix:#filter by columns beginning in prefix
        headings = [i for i in headings if i.startswith(prefix)]
    if experiment:#Experiment name goes on the end
        headings = [i for i in headings if i.endswith(experiment)]
    for req in contains:
        headings = [i for i in headings if req in i]
    for req in not_contains:
        headings = [i for i in headings if not req in i]
    return headings

def load_dataset(file=None, usecols=None, prefix='Abundance:',
                experiment=None, contains=[], not_contains=['corrected', 'count'],
                index='Proteins Unique Sequence ID'):
    #Takes a file and returns a dataframe.
    #    file: the file path to read from
    #    The rest of the paramters are used to select the columns.
    #    By default, it will look for ones starting with 'Reporter intensity'
    #        that do not contain 'count' or 'corrected' and use the 'Protein IDs'
    #        column as the indecies. These will be the raw intensity values.

    #import pdb; pdb.set_trace()
    file = download_file()
    if not file:
        #If no file is named, it will open a filedialog
        file=find_file()
        if not file:#No file selected
            print ("No file selected.")
            return False
    if not usecols:#User has the option of naming columns explicitly=0)
        usecols = get_cols(file, prefix=prefix, experiment=experiment,
                            contains=contains, not_contains=not_contains)
    final_col = [index]
    for i in usecols: final_col.append(i)
    try: df = pd.read_csv(file, sep='\t', header=0, index_col=0, usecols=final_col)
    except: df = pd.read_excel(file, sep='\t', header=0, index_col=0, usecols=final_col)

    return df
