function includeHTML() {
  var z, i, elmnt, file, xhttp;
  z = document.getElementsByTagName("*");
  for (i = 0; i < z.length; i++) {
    elmnt = z[i];
    file = elmnt.getAttribute("w3-include-html");
    if (file) {
      xhttp = new XMLHttpRequest();
      xhttp.onreadystatechange = function() {
        if (this.readyState == 4) {
          if (this.status == 200) {
            elmnt.innerHTML = this.responseText;
            if (elmnt.innerHTML.length === 0) {
              removeAnnouncementBanner();
              return;
            }
          } else if (this.status == 404) {
            removeAnnouncementBanner();
            return;
          }
          elmnt.removeAttribute("w3-include-html");
          includeHTML();
        }
      };
      xhttp.open("GET", file, true);
      xhttp.send();
      return;
    }
  }
};

function removeAnnouncementBanner() {
  var banner = document.getElementsByClassName("announcement");
  banner[0].remove();
};
