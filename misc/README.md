# README for misc

## A home for tools that don't fit in other categories.

---

#### cors.xml

2026.03.31 This was created with help from hyak admin (Nebjosa) to fix a problem getting the LO website to use json particle tracks pulled from kopah. Here is the full solution email:

Based on the error, this is CORS (Cross-Origin Resource Sharing) issue. CORS is a security feature in web browsers which controls how a web page from one origin (domain) can request resources from another origin. To allow your website to access the JSON file, you’ll need to configure a CORS policy on your bucket.

For example, you can create file called cors.xml with the following content:

<CORSConfiguration>
  <CORSRule>
    <AllowedOrigin>https://faculty.washington.edu</AllowedOrigin>
    <AllowedMethod>GET</AllowedMethod>
    <AllowedHeader>*</AllowedHeader>
  </CORSRule>
</CORSConfiguration>

Then, you can apply the CORS policy to your bucket with s3cmd:

s3cmd setcors cors.xml s3://liveocean-web

You can verify with

s3cmd info s3://liveocean-web

Once applied, your JSON file should be accessible from your website without the CORS error.

---

