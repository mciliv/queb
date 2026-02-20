from flask import Flask, render_template, request, jsonify
from astral import LocationInfo
from astral.sun import sun
from datetime import datetime, timezone
import pytz

app = Flask(__name__)

def get_sunset_time(lat, lon):
    city = LocationInfo(latitude=lat, longitude=lon)
    s = sun(city.observer, date=datetime.now().date(), tzinfo=pytz.utc)
    return s['sunset']

@app.route('/')
def index():
    return """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Sunset Tracker</title>
        <script>
            function getLocation() {
                if (navigator.geolocation) {
                    navigator.geolocation.getCurrentPosition(sendPosition);
                } else {
                    alert("Geolocation is not supported by this browser.");
                }
            }
            function sendPosition(position) {
                fetch('/sunset', {
                    method: 'POST',
                    headers: {'Content-Type': 'application/json'},
                    body: JSON.stringify({lat: position.coords.latitude, lon: position.coords.longitude})
                })
                .then(response => response.json())
                .then(data => {
                    document.getElementById("sunsetTime").innerText = "Sunset Time: " + data.sunset_time;
                    document.getElementById("timeUntilSunset").innerText = "Time Until Sunset: " + data.time_until_sunset + " seconds";
                })
                .catch(error => console.error('Error:', error));
            }
        </script>
    </head>
    <body onload="getLocation()">
        <h1>Sunset Tracker</h1>
        <p id="sunsetTime">Fetching sunset time...</p>
        <p id="timeUntilSunset"></p>
    </body>
    </html>
    """

@app.route('/sunset', methods=['POST'])
def sunset():
    data = request.json
    lat = data.get('lat')
    lon = data.get('lon')
    if lat is None or lon is None:
        return jsonify({'error': 'Missing coordinates'}), 400
    
    sunset_time = get_sunset_time(lat, lon)
    sunset_utc = sunset_time.astimezone(timezone.utc)
    now_utc = datetime.now(timezone.utc)
    time_until_sunset = (sunset_utc - now_utc).total_seconds()
    
    return jsonify({
        'sunset_time': sunset_time.strftime('%Y-%m-%d %H:%M:%S %Z'),
        'time_until_sunset': time_until_sunset
    })

if __name__ == '__main__':
    app.run(debug=True)

