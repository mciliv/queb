require 'sinatra'
require 'sinatra/json'
require 'json'
require 'fileutils'
require 'securerandom'

set :port, ENV['PORT'] || 3333
set :public_folder, File.dirname(__FILE__) + '/public'
set :bind, '0.0.0.0'
LOG_FILE = File.join(File.dirname(__FILE__), 'logs', 'access.log')
ADMIN_KEY_FILE = File.join(File.dirname(__FILE__), '.admin_key')

# Ensure logs directory exists
FileUtils.mkdir_p(File.dirname(LOG_FILE))

# Auto-generate admin key if not exists
def get_or_create_admin_key
  return ENV['ADMIN_KEY'] if ENV['ADMIN_KEY']

  if File.exist?(ADMIN_KEY_FILE)
    File.read(ADMIN_KEY_FILE).strip
  else
    new_key = SecureRandom.urlsafe_base64(32)
    File.write(ADMIN_KEY_FILE, new_key)
    new_key
  end
end

ADMIN_KEY = get_or_create_admin_key

# Security headers middleware
before do
  # Prevent clickjacking
  headers 'X-Frame-Options' => 'DENY'

  # Prevent MIME sniffing
  headers 'X-Content-Type-Options' => 'nosniff'

  # XSS protection
  headers 'X-XSS-Protection' => '1; mode=block'

  # Content Security Policy - only allow same origin
  headers 'Content-Security-Policy' => "default-src 'self'; media-src 'self'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline'"

  # Disable caching for sensitive content
  headers 'Cache-Control' => 'no-store, no-cache, must-revalidate, private'
  headers 'Pragma' => 'no-cache'
  headers 'Expires' => '0'
end

# Viewer tracking helper
def track_access(resource_type)
  timestamp = Time.now.utc.iso8601
  ip = request.env['HTTP_X_FORWARDED_FOR'] || request.ip
  user_agent = request.user_agent || 'Unknown'
  referer = request.referer || 'Direct'

  log_entry = {
    timestamp: timestamp,
    resource: resource_type,
    ip: ip,
    userAgent: user_agent,
    referer: referer,
    path: request.path
  }

  log_line = log_entry.to_json + "\n"

  begin
    File.open(LOG_FILE, 'a') do |file|
      file.write(log_line)
    end
  rescue => e
    puts "Failed to write log: #{e.message}"
  end
end

# Upload page
get '/upload' do
  admin_key = params[:key]
  if admin_key != ADMIN_KEY
    status 403
    return 'Access denied. Use ?key=YOUR_ADMIN_KEY'
  end
  send_file File.join(settings.public_folder, 'upload.html')
end

# Handle video upload
post '/upload' do
  admin_key = params[:key]
  if admin_key != ADMIN_KEY
    status 403
    return json({ error: 'Access denied' })
  end

  if params[:video] && params[:video][:tempfile]
    video_path = File.join(settings.public_folder, 'cc472.mp4')
    begin
      File.open(video_path, 'wb') do |f|
        f.write(params[:video][:tempfile].read)
      end
      file_size = File.size(video_path)
      json({ success: true, message: 'Video uploaded successfully', size: file_size })
    rescue => e
      status 500
      json({ error: "Upload failed: #{e.message}" })
    end
  else
    status 400
    json({ error: 'No video file provided' })
  end
end

# Main video page
get '/' do
  track_access('page')
  send_file File.join(settings.public_folder, 'index.html')
end

# Video access tracking
get '/video/cc472.mp4' do
  track_access('video')
  video_path = File.join(settings.public_folder, 'cc472.mp4')

  if File.exist?(video_path)
    send_file video_path, type: 'video/mp4'
  else
    status 404
    'Video not found'
  end
end

# Logs viewer (secured endpoint)
get '/admin/logs' do
  admin_key = params[:key]

  if admin_key != ADMIN_KEY
    status 403
    return 'Access denied. Use ?key=YOUR_ADMIN_KEY'
  end

  send_file File.join(settings.public_folder, 'logs.html')
end

# API endpoint to fetch logs data
get '/admin/api/logs' do
  admin_key = params[:key]

  if admin_key != ADMIN_KEY
    status 403
    return json({ error: 'Access denied' })
  end

  unless File.exist?(LOG_FILE)
    return json([])
  end

  begin
    content = File.read(LOG_FILE)
    logs = content
      .strip
      .split("\n")
      .select { |line| !line.empty? }
      .map do |line|
        begin
          JSON.parse(line)
        rescue JSON::ParserError
          nil
        end
      end
      .compact
      .reverse  # Most recent first

    json logs
  rescue => e
    status 500
    json({ error: 'Failed to read logs' })
  end
end