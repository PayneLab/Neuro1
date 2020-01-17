import os
import requests
import getpass
import bs4

def download_file(download_to_path="data/datafile", url_file_path="data/url.txt", password_file_path="data/password.txt"):
    """Download a file from a given url to the specified location.

    Parameters:
    path (str): The path to the file to save the file to on the local machine.

    Returns:
    str: The path the file was downloaded to.
    """

    import pdb; pdb.set_trace()
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
