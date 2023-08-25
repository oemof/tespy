function removeAnnouncementBanner() {
  var banner = document.getElementsByClassName("announcement")
  banner[0].remove()
};


window.addEventListener('load', async () => {
  const elmnt = document.querySelector('div[oemof-announcement]')
  const response = await fetch(elmnt.getAttribute('oemof-announcement'))
  if (response.status === 200) {
    elmnt.innerHTML = (
      await response.text()
    )
    if (elmnt.innerHTML.length === 0) {
      removeAnnouncementBanner()
    }
  } else if (response.status === 404) {
    removeAnnouncementBanner()
  }
})