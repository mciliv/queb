#!/usr/bin/env node

const puppeteer = require('puppeteer');

async function checkExistingTab() {
  console.log('🔍 Checking for existing Chrome molecular testing tab...\n');

  try {
    // Try to connect to existing Chrome instance
    const browser = await puppeteer.connect({ 
      browserURL: 'http://127.0.0.1:9224',
      defaultViewport: { width: 1600, height: 1000 }
    });
    
    console.log('✅ Found existing Chrome instance!');
    
    const pages = await browser.pages();
    console.log(`📄 Total tabs open: ${pages.length}`);
    
    let molecularTab = null;
    
    for (let i = 0; i < pages.length; i++) {
      try {
        const page = pages[i];
        const url = page.url();
        const title = await page.title();
        
        console.log(`   Tab ${i+1}: ${title} (${url})`);
        
        if (url.includes('localhost:3001') || url.includes('molecular')) {
          molecularTab = page;
          console.log(`   🧪 ← This is your molecular testing tab!`);
        }
        
      } catch (error) {
        console.log(`   Tab ${i+1}: (Unable to access)`);
      }
    }
    
    if (molecularTab) {
      console.log('\n✅ Molecular testing tab found and ready to use!');
      console.log('💡 Run npm run test:persistent-tab to use this existing tab');
    } else {
      console.log('\n📄 No molecular testing tab found');
      console.log('💡 Run npm run test:persistent-tab to create one');
    }
    
    console.log('\n🔧 Chrome debugging info:');
    console.log('   Debug port: 9224');
    console.log('   Debug URL: http://127.0.0.1:9224');
    
    // Don't disconnect - leave browser running
    
  } catch (error) {
    console.log('❌ No existing Chrome instance found');
    console.log(`   Error: ${error.message}`);
    console.log('\n💡 Run npm run test:persistent-tab to create a new instance');
  }
}

checkExistingTab().catch(console.error);